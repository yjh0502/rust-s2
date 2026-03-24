/*
Copyright 2014 Google Inc. All rights reserved.
Copyright 2017 Jhyun Yu. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

use std::cmp::min;

use crate::consts::search_lower_by;
use crate::r3::vector::Vector;
use crate::s1::Angle;
use crate::s2::cap::Cap;
use crate::s2::cell::Cell;
use crate::s2::cellid::*;
use crate::s2::metric::{AVG_AREAMETRIC, MIN_WIDTHMETRIC};
use crate::s2::point::Point;
use crate::s2::rect::Rect;
use crate::s2::region::Region;

/// A CellUnion is a collection of CellIDs.
///
/// It is normalized if it is sorted, and does not contain redundancy.
/// Specifically, it may not contain the same CellID twice, nor a CellID that
/// is contained by another, nor the four sibling CellIDs that are children of
/// a single higher level CellID.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CellUnion(pub Vec<CellID>);

impl CellUnion {
    /// from_range creates a CellUnion that covers the half-open range
    /// of leaf cells [begin, end). If begin == end the resulting union is empty.
    /// This requires that begin and end are both leaves, and begin <= end.
    /// To create a closed-ended range, pass in end.next().
    pub fn from_range(begin: CellID, end: CellID) -> Self {
        let mut v = Vec::new();
        let mut cur = begin.max_tile(&end);
        while cur != end {
            v.push(cur);
            cur = cur.next().max_tile(&end);
        }
        CellUnion(v)
    }

    /// Returns a CellUnion that covers the entire sphere.
    pub fn whole_sphere() -> Self {
        CellUnion(vec![
            CellID::from_face(0),
            CellID::from_face(1),
            CellID::from_face(2),
            CellID::from_face(3),
            CellID::from_face(4),
            CellID::from_face(5),
        ])
    }

    /// normalize normalizes the CellUnion.
    pub fn normalize(&mut self) {
        let v = &mut self.0;
        v.sort();

        let mut output: Vec<CellID> = Vec::with_capacity(v.len());
        for ci in v.iter_mut() {
            // The first two passes here either ignore this new candidate,
            // or remove previously accepted cells that are covered by this candidate.

            // Ignore this cell if it is contained by the previous one.
            // We only need to check the last accepted cell. The ordering of the
            // cells implies containment (but not the converse), and output has no redundancy,
            // so if this candidate is not contained by the last accepted cell
            // then it cannot be contained by any previously accepted cell.
            if let Some(true) = output.last().map(|last| last.contains(ci)) {
                continue;
            }

            // Discard any previously accepted cells contained by this one.
            // This could be any contiguous trailing subsequence, but it can't be
            // a discontiguous subsequence because of the containment property of
            // sorted S2 cells mentioned above.
            while let Some(true) = output.last().map(|last| ci.contains(last)) {
                output.pop();
            }

            // See if the last three cells plus this one can be collapsed.
            // We loop because collapsing three accepted cells and adding a higher level cell
            // could cascade into previously accepted cells.
            let mut ci = *ci;
            while output.len() >= 3 {
                {
                    let fin = &output[(output.len() - 3)..];

                    // fast XOR test; a necessary but not sufficient condition
                    if fin[0].0 ^ fin[1].0 ^ fin[2].0 ^ ci.0 != 0 {
                        break;
                    }

                    // more expensive test; exact.
                    // Compute the two bit mask for the encoded child position,
                    // then see if they all agree.
                    let mut mask = ci.lsb() << 1;
                    mask = !(mask + (mask << 1));
                    let should = ci.0 & mask;
                    if (fin[0].0 & mask != should)
                        || (fin[1].0 & mask != should)
                        || (fin[2].0 & mask != should)
                        || ci.is_face()
                    {
                        break;
                    }
                }

                // output = &output[0..(output.len() - 3)];
                for _ in 0..3 {
                    output.pop();
                }
                ci = ci.immediate_parent();
            }
            output.push(ci);
        }

        // self.0 = output;
        v.clear();
        v.extend_from_slice(&output);
    }

    /// intersects_cellid reports whether this cell union intersects the given cell ID.
    /// This method assumes that the CellUnion has been normalized.
    pub fn intersects_cellid(&self, id: &CellID) -> bool {
        let v = &self.0;
        // Find index of array item that occurs directly after our probe cell:
        let i = search_lower_by(v.len(), |i| id.0 < v[i].0);
        if i != v.len() && v[i].range_min() <= id.range_max() {
            return true;
        }
        i != 0 && v[i - 1].range_max() >= id.range_min()
    }

    /// contains_cellid reports whether the cell union contains the given cell ID.
    /// Containment is defined with respect to regions, e.g. a cell contains its 4 children.
    ///
    /// This method assumes that the CellUnion has been normalized.
    pub fn contains_cellid(&self, id: &CellID) -> bool {
        let v = &self.0;
        // Find index of array item that occurs directly after our probe cell:
        let i = search_lower_by(v.len(), |i| id.0 < v[i].0);
        if i != v.len() && v[i].range_min().0 <= id.0 {
            return true;
        }
        i != 0 && v[i - 1].range_max().0 >= id.0
    }

    /// denormalize replaces this CellUnion with an expanded version of the
    /// CellUnion where any cell whose level is less than minLevel or where
    /// (level - minLevel) is not a multiple of levelMod is replaced by its
    /// children, until either both of these conditions are satisfied or the
    /// maximum level is reached.
    pub fn denormalize(&mut self, min_level: u64, level_mod: u64) {
        let mut v: Vec<CellID> = Vec::new();
        for &id in self.0.iter() {
            let level = id.level();
            let mut new_level = level;
            if new_level < min_level {
                new_level = min_level;
            }
            if level_mod > 1 {
                new_level += (MAX_LEVEL - (new_level - min_level)) % level_mod;
                if new_level > MAX_LEVEL {
                    new_level = MAX_LEVEL;
                }
            }
            if new_level == level {
                v.push(id);
            } else {
                for id in id.child_iter_at_level(new_level) {
                    v.push(id);
                }
            }
        }
        self.0.clear();
        self.0.extend_from_slice(&v);
    }

    /// leaf_cells_covered reports the number of leaf cells covered by this cell union.
    /// This will be no more than 6*2^60 for the whole sphere.
    pub fn leaf_cell_covered(&self) -> u64 {
        let mut num_leaves = 0u64;
        for c in self.0.iter() {
            num_leaves += 1 << ((MAX_LEVEL - c.level()) << 1);
        }
        num_leaves
    }

    /// Reports whether the cell IDs are valid, non-overlapping, and sorted in increasing order.
    pub fn is_valid(&self) -> bool {
        for i in 0..self.0.len() {
            if !self.0[i].is_valid() {
                return false;
            }
            if i > 0 && self.0[i - 1].range_max() >= self.0[i].range_min() {
                return false;
            }
        }
        true
    }

    /// Reports whether the union is valid and no four cells share a common parent.
    pub fn is_normalized(&self) -> bool {
        for i in 0..self.0.len() {
            if !self.0[i].is_valid() {
                return false;
            }
            if i > 0 && self.0[i - 1].range_max() >= self.0[i].range_min() {
                return false;
            }
            if i >= 3 && are_siblings(self.0[i - 3], self.0[i - 2], self.0[i - 1], self.0[i]) {
                return false;
            }
        }
        true
    }

    /// Concatenates several unions and normalizes the result (Go: `CellUnionFromUnion`).
    pub fn merge(unions: &[CellUnion]) -> CellUnion {
        let mut v = Vec::new();
        for u in unions {
            v.extend_from_slice(&u.0);
        }
        let mut cu = CellUnion(v);
        cu.normalize();
        cu
    }

    /// Normalized union of two cell unions.
    pub fn union(a: &Self, b: &Self) -> Self {
        Self::merge(&[a.clone(), b.clone()])
    }

    /// Intersection of two cell unions whose cell ID vectors are sorted in increasing order
    /// (as for a [`CellUnion`] built from normalized input).
    ///
    /// This does not call [`normalize`](Self::normalize) on the result, matching upstream S2:
    /// if both inputs are [`is_normalized`](Self::is_normalized), the result is normalized;
    /// otherwise callers can normalize explicitly if they need a canonical representation.
    ///
    /// In debug builds, asserts that inputs are sorted and that the normalization invariant
    /// above holds (mirroring `ABSL_DCHECK` in the C++ implementation).
    pub fn intersection(a: &Self, b: &Self) -> Self {
        let x = &a.0;
        let y = &b.0;
        debug_assert!(x.is_sorted());
        debug_assert!(y.is_sorted());
        // Mirrors `S2CellUnion::GetIntersection` in google/s2geometry `s2cell_union.cc`.
        let mut cu = Vec::new();
        let mut i = 0usize;
        let mut j = 0usize;
        while i < x.len() && j < y.len() {
            let i_min = x[i].range_min();
            let j_min = y[j].range_min();
            if i_min > j_min {
                if x[i] <= y[j].range_max() {
                    cu.push(x[i]);
                    i += 1;
                } else {
                    // Advance j to the first cell whose range could contain x[i].
                    // lower_bound returns >= j+1 (>= 1), so j-1 is always valid.
                    j = lower_bound_cellids(y, j + 1, y.len(), i_min);
                    if x[i] <= y[j - 1].range_max() {
                        j -= 1;
                    }
                }
            } else if j_min > i_min {
                if y[j] <= x[i].range_max() {
                    cu.push(y[j]);
                    j += 1;
                } else {
                    i = lower_bound_cellids(x, i + 1, x.len(), j_min);
                    if y[j] <= x[i - 1].range_max() {
                        i -= 1;
                    }
                }
            } else if x[i] < y[j] {
                cu.push(x[i]);
                i += 1;
            } else {
                cu.push(y[j]);
                j += 1;
            }
        }
        debug_assert!(cu.is_sorted());
        let out = CellUnion(cu);
        debug_assert!(out.is_normalized() || !a.is_normalized() || !b.is_normalized());
        out
    }

    /// Intersection of a cell union with a single cell ID.
    ///
    /// Does not normalize the result. If `x` is [`is_normalized`](Self::is_normalized), the
    /// result is normalized (same contract as upstream S2). `id` must be valid (debug-asserted).
    pub fn intersection_with_cell_id(x: &Self, id: CellID) -> Self {
        debug_assert!(id.is_valid(), "invalid id: {:?}", id);
        if x.contains_cellid(&id) {
            let out = CellUnion(vec![id]);
            debug_assert!(out.is_normalized() || !x.is_normalized());
            return out;
        }
        let id_max = id.range_max();
        let start = lower_bound_cellids(&x.0, 0, x.0.len(), id.range_min());
        let mut cu = Vec::new();
        let mut k = start;
        while k < x.0.len() && x.0[k] <= id_max {
            cu.push(x.0[k]);
            k += 1;
        }
        let out = CellUnion(cu);
        debug_assert!(out.is_normalized() || !x.is_normalized());
        out
    }

    /// Set difference `x - y`.
    ///
    /// Does not normalize the result. If `x` is [`is_normalized`](Self::is_normalized), the
    /// result is normalized (same contract as upstream S2). Call [`normalize`](Self::normalize)
    /// afterward if you need a canonical union regardless of input shape.
    pub fn difference(x: &Self, y: &Self) -> Self {
        // TODO: This is approximately O(N*log(N)), but could probably
        // use similar techniques as intersection() to be more efficient.
        let mut cu = Vec::new();
        for xid in &x.0 {
            cell_union_difference_internal(&mut cu, *xid, y);
        }
        let out = CellUnion(cu);
        debug_assert!(out.is_normalized() || !x.is_normalized());
        out
    }

    /// Whether this union contains every cell ID of `other` (regions semantics).
    pub fn contains_cell_union(&self, other: &Self) -> bool {
        other.0.iter().all(|id| self.contains_cellid(id))
    }

    /// Whether this union intersects any cell ID of `other`.
    pub fn intersects_cell_union(&self, other: &Self) -> bool {
        self.0.iter().any(|id| other.intersects_cellid(id))
    }

    /// Whether the union contains `p` (point containment uses the containing leaf cell).
    pub fn contains_point(&self, p: &Point) -> bool {
        self.contains_cellid(&CellID::from(p))
    }

    /// Expands the union by adding a rim of cells at `expand_level` around its boundary.
    pub fn expand_at_level(&mut self, expand_level: u64) {
        let level_lsb = lsb_for_level(expand_level);
        let mut output: Vec<CellID> = Vec::new();
        let mut i = self.0.len();
        while i > 0 {
            i -= 1;
            let mut id = self.0[i];
            if id.lsb() < level_lsb {
                id = id.parent(expand_level);
                while i > 0 && id.contains(&self.0[i - 1]) {
                    i -= 1;
                }
            }
            output.push(id);
            output.extend(id.all_neighbors(expand_level));
        }
        self.0 = output;
        self.normalize();
    }

    /// Expands the union so it contains all points within `min_radius` of the region,
    /// using cells at most `max_level_diff` levels above the smallest input cell.
    pub fn expand_by_radius(&mut self, min_radius: Angle, max_level_diff: u64) {
        let mut min_level = MAX_LEVEL;
        for cid in &self.0 {
            min_level = min(min_level, cid.level());
        }
        let radius_level = MIN_WIDTHMETRIC.max_level(min_radius.rad());
        if radius_level == 0 && min_radius.rad() > MIN_WIDTHMETRIC.value(0) {
            self.expand_at_level(0);
        }
        self.expand_at_level(min(min_level + max_level_diff, radius_level));
    }

    /// Average area (accurate within about a factor of 1.7).
    pub fn average_area(&self) -> f64 {
        AVG_AREAMETRIC.value(MAX_LEVEL as u8) * self.leaf_cell_covered() as f64
    }

    /// Approximate area (sum of per-cell approximations).
    pub fn approx_area(&self) -> f64 {
        self.0.iter().map(|id| Cell::from(id).approx_area()).sum()
    }

    /// Exact area (sum of per-cell exact areas).
    pub fn exact_area(&self) -> f64 {
        self.0.iter().map(|id| Cell::from(id).exact_area()).sum()
    }
}

impl Region for CellUnion {
    // cap_bound returns a Cap that bounds this entity.
    fn cap_bound(&self) -> Cap {
        if self.0.is_empty() {
            return Cap::empty();
        }

        // Compute the approximate centroid of the region. This won't produce the
        // bounding cap of minimal area, but it should be close enough.
        let mut centroid = Point(Vector {
            x: 0.,
            y: 0.,
            z: 0.,
        });

        for ci in self.0.iter() {
            let area = AVG_AREAMETRIC.value(ci.level() as u8);
            centroid = centroid + (Point::from(ci) * area);
        }

        if centroid.0.x == 0. && centroid.0.y == 0. && centroid.0.z == 0. {
            centroid = Point::from_coords(1., 0., 0.);
        } else {
            centroid = Point(centroid.0.normalize());
        }

        // Use the centroid as the cap axis, and expand the cap angle so that it
        // contains the bounding caps of all the individual cells.  Note that it is
        // *not* sufficient to just bound all the cell vertices because the bounding
        // cap may be concave (i.e. cover more than one hemisphere).
        let mut cap = Cap::from(&centroid);
        for ci in self.0.iter() {
            cap = cap + Cell::from(ci).cap_bound();
        }
        cap
    }

    /// rect_bound returns a Rect that bounds this entity.
    fn rect_bound(&self) -> Rect {
        let mut bound = Rect::empty();
        for c in self.0.iter() {
            bound = bound.union(&Cell::from(c).rect_bound());
        }
        bound
    }

    // contains_cell reports whether this cell union contains the given cell.
    fn contains_cell(&self, c: &Cell) -> bool {
        self.contains_cellid(&c.id)
    }

    // intersects_cell reports whether this cell union intersects the given cell.
    fn intersects_cell(&self, c: &Cell) -> bool {
        self.intersects_cellid(&c.id)
    }

    fn cell_union_bound(&self) -> Vec<CellID> {
        self.cap_bound().cell_union_bound()
    }
}

fn lower_bound_cellids(cells: &[CellID], begin: usize, end: usize, id: CellID) -> usize {
    for i in begin..end {
        if cells[i] >= id {
            return i;
        }
    }
    end
}

fn are_siblings(a: CellID, b: CellID, c: CellID, d: CellID) -> bool {
    // A necessary (but not sufficient) condition is that the XOR of the
    // four cells must be zero.  This is also very fast to test.
    if a.0 ^ b.0 ^ c.0 != d.0 {
        return false;
    }

    // Now we do a slightly more expensive but exact test.  First, compute a
    // mask that blocks out the two bits that encode the child position of
    // `id` with respect to its parent, then check that the other three
    // children all agree with `mask`.
    let mut mask = d.lsb() << 1;
    mask = !(mask + (mask << 1));
    let id_masked = d.0 & mask;
    (a.0 & mask) == id_masked
        && (b.0 & mask) == id_masked
        && (c.0 & mask) == id_masked
        && !d.is_face()
}

fn cell_union_difference_internal(out: &mut Vec<CellID>, id: CellID, other: &CellUnion) {
    // Mirrors `GetDifferenceInternal` in google/s2geometry `s2cell_union.cc`.
    // Add the difference between `id` and `other` to `out`.
    // If they intersect but the difference is non-empty, divide and conquer.
    if !other.intersects_cellid(&id) {
        out.push(id);
        return;
    }
    if !other.contains_cellid(&id) {
        let mut child = id.child_begin();
        let mut i = 0;
        loop {
            cell_union_difference_internal(out, child, other);
            if i == 3 {
                break;
            }
            child = child.next();
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::EPSILON;

    #[test]
    fn test_cellunion_normalization() {
        let mut cu = CellUnion(vec![
            CellID(0x80855c0000000000), // A: a cell over Pittsburg CA
            CellID(0x80855d0000000000), // B, a child of A
            CellID(0x8085634000000000), // first child of X, disjoint from A
            CellID(0x808563c000000000), // second child of X
            CellID(0x80855dc000000000), // a child of B
            CellID(0x808562c000000000), // third child of X
            CellID(0x8085624000000000), // fourth child of X
            CellID(0x80855d0000000000), // B again
        ]);

        let exp = CellUnion(vec![
            CellID(0x80855c0000000000), // A
            CellID(0x8085630000000000), // X
        ]);

        cu.normalize();

        assert_eq!(cu, exp);

        cu.0.push(CellID(0x808562c000000000));
        cu.normalize();
        assert_eq!(cu, exp);
    }

    #[test]
    fn test_cellunion_basic() {
        let mut empty = CellUnion(vec![]);
        empty.normalize();
        assert_eq!(empty.0.len(), 0);

        let face1_id = CellID::from_face(1);
        let face1_cell = Cell::from(&face1_id);
        let mut face1_union = CellUnion(vec![face1_id.clone()]);
        face1_union.normalize();
        assert_eq!(face1_union.0.len(), 1);
        assert_eq!(face1_id, face1_union.0[0]);
        assert_eq!(true, face1_union.contains_cell(&face1_cell));

        let face2_id = CellID::from_face(2);
        let face2_cell = Cell::from(&face2_id);
        let mut face2_union = CellUnion(vec![face2_id.clone()]);
        face2_union.normalize();
        assert_eq!(face2_union.0.len(), 1);
        assert_eq!(face2_id, face2_union.0[0]);
        assert_eq!(true, face2_union.contains_cell(&face2_cell));

        assert_eq!(false, face1_union.contains_cell(&face2_cell));
    }

    fn test_cellunion_case(
        cells: &[CellID],
        contained: &[CellID],
        overlaps: &[CellID],
        disjoint: &[CellID],
    ) {
        let mut v = Vec::with_capacity(cells.len());
        v.extend_from_slice(cells);
        let mut union = CellUnion(v);
        union.normalize();

        // Ensure self-containment tests are correct.
        for id in cells {
            assert_eq!(true, union.intersects_cellid(&id));
            assert_eq!(true, union.contains_cellid(&id));
        }

        // Test for containment specified in test case.
        for id in contained {
            assert_eq!(true, union.intersects_cellid(&id));
            assert_eq!(true, union.contains_cellid(&id));
        }

        // Make sure the CellUnion intersect these cells but do not contain.
        for id in overlaps {
            assert_eq!(true, union.intersects_cellid(&id));
            assert_eq!(false, union.contains_cellid(&id));
        }

        // Negative cases make sure the CellUnion neither contain nor intersect these cells
        for id in disjoint {
            assert_eq!(false, union.intersects_cellid(&id));
            assert_eq!(false, union.contains_cellid(&id));
        }
    }

    #[test]
    fn test_cellunion() {
        test_cellunion_case(
            &[
                // Single cell around NYC, and some simple nearby probes
                CellID(0x89c25c0000000000),
            ],
            &[
                CellID(0x89c25c0000000000).child_begin(),
                CellID(0x89c25c0000000000).child_begin_at_level(28),
            ],
            &[
                CellID(0x89c25c0000000000).immediate_parent(),
                // the whole face
                CellID::from_face(CellID(0x89c25c0000000000).face() as u64),
            ],
            &[
                // Cell next to this one at same level
                CellID(0x89c25c0000000000).next(),
                // Cell next to this one at deep level
                CellID(0x89c25c0000000000).next().child_begin_at_level(28),
                // Big(er) neighbor cell
                CellID(0x89c2700000000000),
                // Very big next door cell.
                CellID(0x89e9000000000000),
                // Very big cell, smaller value than probe
                CellID(0x89c1000000000000),
            ],
        );

        test_cellunion_case(
            &[
                // NYC and SFO:
                CellID(0x89c25b0000000000), // NYC
                CellID(0x89c2590000000000), // NYC
                CellID(0x89c2f70000000000), // NYC
                CellID(0x89c2f50000000000), // NYC
                CellID(0x8085870000000000), // SFO
                CellID(0x8085810000000000), // SFO
                CellID(0x808f7d0000000000), // SFO
                CellID(0x808f7f0000000000), // SFO
            ],
            &[
                CellID(0x808f7ef300000000), // SFO
                CellID(0x808f7e5cf0000000), // SFO
                CellID(0x808587f000000000), // SFO
                CellID(0x89c25ac000000000), // NYC
                CellID(0x89c259a400000000), // NYC
                CellID(0x89c258fa10000000), // NYC
                CellID(0x89c258f174007000), // NYC
            ],
            &[
                CellID(0x808c000000000000), // Big SFO
                CellID(0x89c4000000000000), // Big NYC
            ],
            &[
                CellID(0x89c15a4fcb1bb000), // outside NYC
                CellID(0x89c15a4e4aa95000), // outside NYC
                CellID(0x8094000000000000), // outside SFO (big)
                CellID(0x8096f10000000000), // outside SFO (smaller)
                CellID(0x87c0000000000000), // Midwest very big
            ],
        );

        test_cellunion_case(
            &[
                CellID(0x8100000000000000), // starting around california
                CellID(0x8740000000000000), // adjacent cells at increasing
                CellID(0x8790000000000000), // levels, moving eastward.
                CellID(0x87f4000000000000),
                CellID(0x87f9000000000000), // going down across the midwest
                CellID(0x87ff400000000000),
                CellID(0x87ff900000000000),
                CellID(0x87fff40000000000),
                CellID(0x87fff90000000000),
                CellID(0x87ffff4000000000),
                CellID(0x87ffff9000000000),
                CellID(0x87fffff400000000),
                CellID(0x87fffff900000000),
                CellID(0x87ffffff40000000),
                CellID(0x87ffffff90000000),
                CellID(0x87fffffff4000000),
                CellID(0x87fffffff9000000),
                CellID(0x87ffffffff400000), // to a very small cell in Wisconsin
            ],
            &[
                CellID(0x808f400000000000),
                CellID(0x80eb118b00000000),
                CellID(0x8136a7a11d000000),
                CellID(0x8136a7a11dac0000),
                CellID(0x876c7c0000000000),
                CellID(0x87f96d0000000000),
                CellID(0x87ffffffff400000),
            ],
            &[
                CellID(0x8100000000000000).immediate_parent(),
                CellID(0x8740000000000000).immediate_parent(),
            ],
            &[
                CellID(0x52aaaaaaab300000),
                CellID(0x52aaaaaaacd00000),
                CellID(0x87fffffffa100000),
                CellID(0x87ffffffed500000),
                CellID(0x87ffffffa0100000),
                CellID(0x87fffffed5540000),
                CellID(0x87fffffed6240000),
                CellID(0x52aaaacccb340000),
                CellID(0x87a0000400000000),
                CellID(0x87a000001f000000),
                CellID(0x87a0000029d00000),
                CellID(0x9500000000000000),
            ],
        );
    }

    #[test]
    fn test_cellunion_merge_intersection_difference() {
        let a = CellUnion(vec![CellID::from_face(0).child_begin_at_level(2)]);
        let b = CellUnion(vec![CellID::from_face(1).child_begin_at_level(2)]);
        let m = CellUnion::merge(&[a.clone(), b.clone()]);
        assert_eq!(m.0.len(), 2);
        assert!(m.is_normalized());

        let ab = CellUnion::intersection(&a, &b);
        assert!(ab.0.is_empty());

        let p = CellID::from_face(0).child_begin_at_level(2);
        let ch = p.children();
        let x = CellUnion(vec![p]);
        let y = CellUnion(vec![ch[0]]);
        let d = CellUnion::difference(&x, &y);
        assert_eq!(d.0.len(), 3);
        let merged_back = CellUnion::merge(&[d, y]);
        assert_eq!(merged_back.0, vec![p]);
    }

    #[test]
    fn test_cellunion_contains_intersects_union() {
        let p = CellID::from_face(0).child_begin_at_level(1);
        let ch = p.children();
        let children_u = CellUnion(vec![ch[0], ch[1], ch[2], ch[3]]);
        let parent_u = CellUnion(vec![p]);
        assert!(parent_u.contains_cell_union(&children_u));
        assert!(children_u.intersects_cell_union(&parent_u));
        assert!(!children_u.contains_cell_union(&parent_u));
    }

    #[test]
    fn test_cellunion_intersection_with_cell_id() {
        let p = CellID::from_face(3).child_begin_at_level(4);
        let inside = p.child_begin_at_level(8);
        let mut u = CellUnion(vec![p]);
        u.normalize();
        let cap = CellUnion::intersection_with_cell_id(&u, inside);
        assert_eq!(cap.0, vec![inside]);
    }

    #[test]
    fn test_cellunion_intersection_with_cell_id_parent_covers_multiple_cells() {
        // `x` holds two children of `parent`; `x` does not contain the parent id as one cell,
        // so we exercise the range scan path (not the `contains_cellid` fast path).
        let parent = CellID::from_face(2).child_begin_at_level(5);
        let ch = parent.children();
        let mut x = CellUnion(vec![ch[0], ch[1]]);
        x.normalize();
        assert!(x.is_normalized());

        let cap = CellUnion::intersection_with_cell_id(&x, parent);
        assert_eq!(cap.0, vec![ch[0], ch[1]]);
        assert!(cap.is_normalized());
    }

    #[test]
    fn test_normalized_intersection_and_difference_remain_normalized() {
        let a = CellID::from_face(0).child_begin_at_level(3);
        let b = CellID::from_face(1).child_begin_at_level(3);
        let mut ua = CellUnion(vec![a]);
        let mut ub = CellUnion(vec![b]);
        ua.normalize();
        ub.normalize();

        let inter = CellUnion::intersection(&ua, &ub);
        assert!(inter.is_normalized());
        assert!(inter.0.is_empty());

        let p = CellID::from_face(4).child_begin_at_level(4);
        let ch = p.children();
        let mut ux = CellUnion(vec![p]);
        let mut uy = CellUnion(vec![ch[0]]);
        ux.normalize();
        uy.normalize();
        let diff = CellUnion::difference(&ux, &uy);
        assert!(diff.is_normalized());
        assert_eq!(diff.0.len(), 3);
    }

    #[test]
    fn test_cellunion_from_range_empties_and_full_sphere() {
        let id_begin = CellID::from_face(0).child_begin_at_level(MAX_LEVEL);
        assert!(CellUnion::from_range(id_begin, id_begin).0.is_empty());

        let id_end = CellID::from_face(5).child_end_at_level(MAX_LEVEL);
        assert!(CellUnion::from_range(id_end, id_end).0.is_empty());

        let cu = CellUnion::from_range(id_begin, id_end);
        assert_eq!(cu.0.len(), 6);
        assert!(cu.0.iter().all(|c| c.is_face()));
    }

    #[test]
    fn test_cellunion_leaf_cells_covered_disjoint_mix() {
        let have = vec![
            CellID::from_face(0).child_begin_at_level(MAX_LEVEL),
            CellID::from_face(0),
            CellID::from_face(1).child_begin_at_level(1),
            CellID::from_face(2).child_begin_at_level(2),
            CellID::from_face(2).child_end_at_level(2).prev(),
            CellID::from_face(3).child_begin_at_level(14),
            CellID::from_face(4).child_begin_at_level(27),
            CellID::from_face(4).child_end_at_level(15).prev(),
            CellID::from_face(5).child_begin_at_level(30),
        ];
        let mut cu = CellUnion(have);
        cu.normalize();
        let want = 1u64 + (1 << 6) + (1 << 30) + (1 << 32) + (2 << 56) + (1 << 58) + (1 << 60);
        assert_eq!(cu.leaf_cell_covered(), want);
    }

    #[test]
    fn test_cellunion_contains_point_and_areas() {
        use crate::s1::{Angle, Rad};
        use crate::s2::metric::AVG_AREAMETRIC;

        let id = CellID::from_face(2).child_begin_at_level(10);
        let p = Cell::from(&id).center();
        let cu = CellUnion(vec![id]);
        assert!(cu.contains_point(&p));

        let cell = Cell::from(&id);
        assert_f64_eq!(cu.approx_area(), cell.approx_area());
        assert_f64_eq!(cu.exact_area(), cell.exact_area());
        assert_f64_eq!(
            cu.average_area(),
            AVG_AREAMETRIC.value(MAX_LEVEL as u8) * cu.leaf_cell_covered() as f64
        );

        let mut expanded = cu.clone();
        expanded.expand_by_radius(Angle::from(Rad(1e-6)), 2);
        assert!(expanded.contains_point(&p));
    }

    #[test]
    fn test_cellunion_is_valid_and_is_normalized() {
        // A properly normalized union is both valid and normalized.
        let mut cu = CellUnion(vec![
            CellID::from_face(0).child_begin_at_level(2),
            CellID::from_face(1).child_begin_at_level(3),
        ]);
        cu.normalize();
        assert!(cu.is_valid());
        assert!(cu.is_normalized());

        // Unsorted cells: not valid.
        let unsorted = CellUnion(vec![
            CellID::from_face(2).child_begin_at_level(2),
            CellID::from_face(0).child_begin_at_level(2),
        ]);
        assert!(!unsorted.is_valid());

        // Overlapping cells (parent + child): not valid.
        let parent = CellID::from_face(0).child_begin_at_level(2);
        let child = parent.child_begin_at_level(4);
        let overlapping = CellUnion(vec![parent, child]);
        assert!(!overlapping.is_valid());

        // Four siblings without collapsing: valid but not normalized.
        let p = CellID::from_face(3).child_begin_at_level(5);
        let ch = p.children();
        let siblings = CellUnion(vec![ch[0], ch[1], ch[2], ch[3]]);
        assert!(siblings.is_valid());
        assert!(!siblings.is_normalized());

        // Empty union is both valid and normalized.
        let empty = CellUnion(vec![]);
        assert!(empty.is_valid());
        assert!(empty.is_normalized());
    }

    #[test]
    fn test_cellunion_intersection_overlapping() {
        // Parent in one union, children in the other — real overlap.
        let parent = CellID::from_face(0).child_begin_at_level(3);
        let ch = parent.children();

        let u_parent = CellUnion(vec![parent]);
        // Four siblings are valid but not normalized; C++ likewise allows
        // GetIntersection without normalizing the output when an input is not.
        let u_children = CellUnion(vec![ch[0], ch[1], ch[2], ch[3]]);

        // Intersection yields the four leaf cells (same as upstream before any
        // second Normalize pass); collapsing to the parent requires normalize().
        let inter = CellUnion::intersection(&u_parent, &u_children);
        assert_eq!(inter.0.len(), 4);
        for c in &ch {
            assert!(inter.contains_cellid(c));
        }
        let mut inter_norm = inter.clone();
        inter_norm.normalize();
        assert_eq!(inter_norm.0.len(), 1);
        assert_eq!(inter_norm.0[0], parent);

        // Intersect with just two children.
        let u_two = CellUnion(vec![ch[0], ch[2]]);
        let inter2 = CellUnion::intersection(&u_parent, &u_two);
        assert_eq!(inter2.0.len(), 2);
        assert!(inter2.contains_cellid(&ch[0]));
        assert!(inter2.contains_cellid(&ch[2]));
        assert!(!inter2.contains_cellid(&ch[1]));

        // Intersection is commutative.
        let inter3 = CellUnion::intersection(&u_two, &u_parent);
        assert_eq!(inter2, inter3);
    }

    #[test]
    fn test_cellunion_expand_at_level() {
        // A single cell expanded at its own level should include itself + its neighbors.
        let id = CellID::from_face(1).child_begin_at_level(5);
        let mut cu = CellUnion(vec![id]);
        cu.normalize();
        cu.expand_at_level(5);

        // Must still contain the original cell.
        assert!(cu.contains_cellid(&id));

        // Must contain all edge neighbors at the same level.
        for nbr in &id.edge_neighbors() {
            assert!(
                cu.contains_cellid(nbr),
                "expanded union should contain edge neighbor {}",
                nbr
            );
        }

        // Expanding a deeper cell at a coarser level: cell gets promoted to its parent
        // at that level, so we also get the parent's neighbors.
        let deep = CellID::from_face(2).child_begin_at_level(10);
        let mut cu2 = CellUnion(vec![deep]);
        cu2.normalize();
        cu2.expand_at_level(8);
        let promoted = deep.parent(8);
        assert!(cu2.contains_cellid(&promoted));
        for nbr in &promoted.edge_neighbors() {
            assert!(cu2.contains_cellid(nbr));
        }
    }

    #[test]
    fn test_cellunion_difference_partial_overlap() {
        let parent = CellID::from_face(0).child_begin_at_level(4);
        let ch = parent.children();
        // x has child 0 and child 1
        let x = CellUnion(vec![ch[0], ch[1]]);
        // y has child 1 and child 2
        let y = CellUnion(vec![ch[1], ch[2]]);

        // x - y should be just child 0
        let diff = CellUnion::difference(&x, &y);
        assert_eq!(diff.0, vec![ch[0]]);

        // y - x should be just child 2
        let diff2 = CellUnion::difference(&y, &x);
        assert_eq!(diff2.0, vec![ch[2]]);

        // (x - y) ∪ (x ∩ y) should reconstruct x
        let inter = CellUnion::intersection(&x, &y);
        assert_eq!(inter.0, vec![ch[1]]);
        let reconstructed = CellUnion::merge(&[diff, inter]);
        assert_eq!(reconstructed.0, vec![ch[0], ch[1]]);
    }

    #[test]
    fn test_cellunion_whole_sphere() {
        let cu = CellUnion::whole_sphere();
        assert_eq!(cu.0.len(), 6);
        assert!(cu.is_valid());
        assert!(cu.is_normalized());
        for face in 0..6 {
            let id = CellID::from_face(face);
            assert!(cu.contains_cellid(&id));
        }
        // Must contain an arbitrary deeply nested cell.
        let deep = CellID::from_face(3).child_begin_at_level(MAX_LEVEL);
        assert!(cu.contains_cellid(&deep));
    }

    #[test]
    fn test_cellunion_empty() {
        let empty = CellUnion(vec![]);

        // Valid and normalized.
        assert!(empty.is_valid());
        assert!(empty.is_normalized());

        // Contains / intersects nothing.
        let face0 = CellID::from_face(0);
        assert!(!empty.contains_cellid(&face0));
        assert!(!empty.intersects_cellid(&face0));

        let p = Point::from(&face0);
        assert!(!empty.contains_point(&p));

        // Operations with empty unions.
        let other = CellUnion(vec![face0]);
        assert!(!empty.contains_cell_union(&other));
        assert!(!empty.intersects_cell_union(&other));
        assert!(empty.contains_cell_union(&CellUnion(vec![])));

        // Union / intersection / difference with empty.
        let u = CellUnion::union(&empty, &other);
        assert_eq!(u.0, vec![face0]);
        let u2 = CellUnion::union(&other, &empty);
        assert_eq!(u2.0, vec![face0]);

        let inter = CellUnion::intersection(&empty, &other);
        assert!(inter.0.is_empty());

        let diff = CellUnion::difference(&other, &empty);
        assert_eq!(diff.0, vec![face0]);
        let diff2 = CellUnion::difference(&empty, &other);
        assert!(diff2.0.is_empty());

        // Areas.
        assert_eq!(empty.average_area(), 0.0);
        assert_eq!(empty.approx_area(), 0.0);
        assert_eq!(empty.exact_area(), 0.0);

        // Leaf cells covered.
        assert_eq!(empty.leaf_cell_covered(), 0);

        // Region bounds.
        let cap = empty.cap_bound();
        assert!(cap.is_empty());
    }

    #[test]
    fn test_cellunion_invalid_cellid_not_valid() {
        // A CellUnion containing an invalid CellID(0) should not be valid.
        let cu = CellUnion(vec![CellID(0)]);
        assert!(!cu.is_valid());
        assert!(!cu.is_normalized());
    }

    #[test]
    fn test_cellunion_duplicate_cells_not_valid() {
        let id = CellID::from_face(1).child_begin_at_level(5);
        let cu = CellUnion(vec![id, id]);
        assert!(!cu.is_valid());
    }

    #[test]
    fn test_are_siblings() {
        // The four children of a cell are siblings.
        let parent = CellID::from_face(2).child_begin_at_level(4);
        let ch = parent.children();
        assert!(are_siblings(ch[0], ch[1], ch[2], ch[3]));

        // Permutation should still work since the XOR test is order-sensitive
        // but the mask test checks identity regardless of order.
        assert!(are_siblings(ch[3], ch[2], ch[1], ch[0]));

        // Replacing one child with an unrelated cell makes them non-siblings.
        let other = CellID::from_face(3).child_begin_at_level(5);
        assert!(!are_siblings(ch[0], ch[1], ch[2], other));

        // Face cells are never siblings (guard in are_siblings).
        let faces: Vec<CellID> = (0..4).map(CellID::from_face).collect();
        assert!(!are_siblings(faces[0], faces[1], faces[2], faces[3]));
    }

    #[test]
    fn test_cellunion_leaf_cells_covered_full_sphere() {
        // 5-face union: all faces except one.
        let five_faces = CellUnion((0..5).map(CellID::from_face).collect());
        let one_face_leaves = 1u64 << (MAX_LEVEL * 2);
        assert_eq!(five_faces.leaf_cell_covered(), 5 * one_face_leaves);

        // Full sphere (6 faces).
        let full = CellUnion::whole_sphere();
        assert_eq!(full.leaf_cell_covered(), 6 * one_face_leaves);

        // Expanding a single face at level 0 adds its 4 edge-neighbor faces.
        // Face 0's opposite face (face 5) shares no edges, so we get 5 faces.
        let id = CellID::from_face(0);
        let mut expanded = CellUnion(vec![id]);
        expanded.expand_at_level(0);
        assert_eq!(expanded.leaf_cell_covered(), 5 * one_face_leaves);
    }

    fn test_denorm_case(
        name: &str,
        min_level: u64,
        level_mod: u64,
        mut cu: CellUnion,
        exp: CellUnion,
    ) {
        cu.denormalize(min_level, level_mod);
        assert_eq!(exp, cu, "{}", name);
    }

    fn cfbl(face: u64, level: u64) -> CellID {
        CellID::from_face(face).child_begin_at_level(level)
    }

    #[test]
    fn test_cellunion_denormalize() {
        test_denorm_case(
            "not expanded, level mod == 1",
            10,
            1,
            CellUnion(vec![cfbl(2, 11), cfbl(2, 11), cfbl(3, 14), cfbl(0, 10)]),
            CellUnion(vec![cfbl(2, 11), cfbl(2, 11), cfbl(3, 14), cfbl(0, 10)]),
        );

        test_denorm_case(
            "not expanded, level mod > 1",
            10,
            2,
            CellUnion(vec![cfbl(2, 12), cfbl(2, 12), cfbl(3, 14), cfbl(0, 10)]),
            CellUnion(vec![cfbl(2, 12), cfbl(2, 12), cfbl(3, 14), cfbl(0, 10)]),
        );

        test_denorm_case(
            "expended, (level - min_level) is not multiple of level mod",
            10,
            3,
            CellUnion(vec![cfbl(2, 12), cfbl(5, 11)]),
            CellUnion(vec![
                cfbl(2, 12).children()[0],
                cfbl(2, 12).children()[1],
                cfbl(2, 12).children()[2],
                cfbl(2, 12).children()[3],
                cfbl(5, 11).children()[0].children()[0],
                cfbl(5, 11).children()[0].children()[1],
                cfbl(5, 11).children()[0].children()[2],
                cfbl(5, 11).children()[0].children()[3],
                cfbl(5, 11).children()[1].children()[0],
                cfbl(5, 11).children()[1].children()[1],
                cfbl(5, 11).children()[1].children()[2],
                cfbl(5, 11).children()[1].children()[3],
                cfbl(5, 11).children()[2].children()[0],
                cfbl(5, 11).children()[2].children()[1],
                cfbl(5, 11).children()[2].children()[2],
                cfbl(5, 11).children()[2].children()[3],
                cfbl(5, 11).children()[3].children()[0],
                cfbl(5, 11).children()[3].children()[1],
                cfbl(5, 11).children()[3].children()[2],
                cfbl(5, 11).children()[3].children()[3],
            ]),
        );

        test_denorm_case(
            "expended, level < min_level",
            10,
            3,
            CellUnion(vec![cfbl(2, 9)]),
            CellUnion(vec![
                cfbl(2, 9).children()[0],
                cfbl(2, 9).children()[1],
                cfbl(2, 9).children()[2],
                cfbl(2, 9).children()[3],
            ]),
        );
    }
}

// Binary Encode/Decode for CellUnion (as in golang/geo) is not implemented.
