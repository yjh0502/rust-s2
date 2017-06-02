/*
Copyright 2014 Google Inc. All rights reserved.
Copyright 2017 Jihyun Yu. All rights reserved.

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


use std;
use std::u64;
use std::iter::*;

use consts::clamp;
use r1;
use r2;
use r3::vector::Vector;
use s1::angle::Angle;
use s2::stuv::*;
use s2::point::Point;
use s2::latlng::*;

/// CellID uniquely identifies a cell in the S2 cell decomposition.
/// The most significant 3 bits encode the face number (0-5). The
/// remaining 61 bits encode the position of the center of this cell
/// along the Hilbert curve on that face. The zero value and the value
/// (1<<64)-1 are invalid cell IDs. The first compares less than any
/// valid cell ID, the second as greater than any valid cell ID.
///
/// Sequentially increasing cell IDs follow a continuous space-filling curve
/// over the entire sphere. They have the following properties:
///
///  - The ID of a cell at level k consists of a 3-bit face number followed
///    by k bit pairs that recursively select one of the four children of
///    each cell. The next bit is always 1, and all other bits are 0.
///    Therefore, the level of a cell is determined by the position of its
///    lowest-numbered bit that is turned on (for a cell at level k, this
///    position is 2 * (maxLevel - k)).
///
///  - The ID of a parent cell is at the midpoint of the range of IDs spanned
///    by its children (or by its descendants at any level).
///
/// Leaf cells are often used to represent points on the unit sphere, and
/// this type provides methods for converting directly between these two
/// representations. For cells that represent 2D regions rather than
/// discrete point, it is better to use Cells.
#[derive(Clone,Copy,PartialEq,Eq,PartialOrd,Ord,Hash)]
pub struct CellID(pub u64);

const FACE_BITS: u64 = 3;
pub const NUM_FACES: u8 = 6;
pub const MAX_LEVEL: u64 = 30;
pub const POS_BITS: u64 = (2 * MAX_LEVEL) + 1;
pub const MAX_SIZE: u64 = 1 << MAX_LEVEL;
const WRAP_OFFSET: u64 = (NUM_FACES as u64) << POS_BITS;

const MAX_SIZE_I32: i32 = MAX_SIZE as i32;
const MAX_SIZE_F64: f64 = MAX_SIZE as f64;

const LOOKUP_BITS: u64 = 4;
const INVERT_MASK: u8 = 0x02;

fn lsb_for_level(level: u64) -> u64 {
    1 << (2 * (MAX_LEVEL - level))
}

//TODO private
pub fn size_ij(level: u64) -> u64 {
    1 << (MAX_LEVEL - level)
}

impl CellID {
    /// from_pos_level returns a cell given its face in the range
    /// [0,5], the 61-bit Hilbert curve position pos within that face, and
    /// the level in the range [0,maxLevel]. The position in the cell ID
    /// will be truncated to correspond to the Hilbert curve position at
    /// the center of the returned cell.
    pub fn from_face_pos_level(face: u64, pos: u64, level: u64) -> Self {
        CellID((face << POS_BITS) + (pos | 1)).parent(level)
    }

    /// from_face returns the cell corresponding to a given S2 cube face.
    pub fn from_face(face: u64) -> Self {
        CellID((face << POS_BITS) + lsb_for_level(0))
    }

    /// from_ij returns a leaf cell given its cube face (range 0..5) and IJ coordinates.
    fn from_face_ij_wrap(face: u8, mut i: i32, mut j: i32) -> Self {
        // Convert i and j to the coordinates of a leaf cell just beyond the
        // boundary of this face.  This prevents 32-bit overflow in the case
        // of finding the neighbors of a face cell.
        i = clamp(i, -1i32, MAX_SIZE_I32);
        j = clamp(j, -1i32, MAX_SIZE_I32);

        const SCALE: f64 = 1.0 / (MAX_SIZE as f64);
        const LIMIT: f64 = 1f64 + std::f64::EPSILON;

        let u = clamp(SCALE * (2. * (i as f64) + 1. - MAX_SIZE_F64), -LIMIT, LIMIT);
        let v = clamp(SCALE * (2. * (j as f64) + 1. - MAX_SIZE_F64), -LIMIT, LIMIT);

        // Find the leaf cell coordinates on the adjacent face, and convert
        // them to a cell id at the appropriate level.
        let (f, u, v) = xyz_to_face_uv(&face_uv_to_xyz(face, u, v));
        Self::from_face_ij(f, st_to_ij(0.5 * (u + 1.)), st_to_ij(0.5 * (v + 1.)))
    }

    //TODO private
    pub fn from_face_ij(f: u8, i: i32, j: i32) -> Self {
        let mut n = (f as u64) << (POS_BITS - 1);
        let mut bits = (f & SWAP_MASK) as i32;

        let mut k = 7;
        let mask = (1 << LOOKUP_BITS) - 1;
        loop {
            bits += ((i >> (k * LOOKUP_BITS)) & mask) << (LOOKUP_BITS + 2);
            bits += ((j >> (k * LOOKUP_BITS)) & mask) << 2;
            bits = LOOKUP_POS[bits as usize] as i32;
            n |= ((bits >> 2) as u64) << ((k * 2 * LOOKUP_BITS) as u64);
            bits &= (SWAP_MASK | INVERT_MASK) as i32;

            if k == 0 {
                break;
            }
            k -= 1;
        }
        CellID(n * 2 + 1)
    }

    fn from_face_ij_same(f: u8, i: i32, j: i32, same_face: bool) -> Self {
        if same_face {
            Self::from_face_ij(f, i, j)
        } else {
            Self::from_face_ij_wrap(f, i, j)
        }
    }


    /// from_token returns a cell given a hex-encoded string of its uint64 ID.
    pub fn from_token(s: &str) -> CellID {
        match u64::from_str_radix(s, 16) {
            Err(_) => CellID(0),
            Ok(mut v) => {
                if s.len() < 16 {
                    v = v << (4 * (16 - s.len()));
                }
                CellID(v)
            }
        }
    }

    /// to_token returns a hex-encoded string of the uint64 cell id, with leading
    /// zeros included but trailing zeros stripped.
    pub fn to_token(&self) -> String {
        if self.0 == 0 {
            "X".into()
        } else {
            format!("{:016x}", self.0).trim_right_matches('0').into()
        }
    }

    /// is_valid reports whether ci represents a valid cell.
    pub fn is_valid(&self) -> bool {
        self.face() < NUM_FACES && (self.lsb() & 0x1555555555555555 != 0)
    }

    /// face returns the cube face for this cell ID, in the range [0,5].
    pub fn face(&self) -> u8 {
        (self.0 >> POS_BITS) as u8
    }

    /// pos returns the position along the Hilbert curve of this cell ID, in the range [0,2^posBits-1].
    pub fn pos(&self) -> u64 {
        self.0 & ((!0u64) >> FACE_BITS)
    }

    /// level returns the subdivision level of this cell ID, in the range [0, maxLevel].
    pub fn level(&self) -> u64 {
        MAX_LEVEL - (self.0.trailing_zeros() >> 1) as u64
    }

    /// is_leaf returns whether this cell ID is at the deepest level;
    /// that is, the level at which the cells are smallest.
    pub fn is_leaf(&self) -> bool {
        self.0 & 1 != 0
    }

    /// child_position returns the child position (0..3) of this cell's
    /// ancestor at the given level, relative to its parent.  The argument
    /// should be in the range 1..kMaxLevel.  For example,
    /// ChildPosition(1) returns the position of this cell's level-1
    /// ancestor within its top-level face cell.
    pub fn child_position(&self, level: u64) -> u64 {
        (self.0 >> (2 * (MAX_LEVEL - level) + 1)) & 3
    }

    /// parent returns the cell at the given level, which must be no greater than the current level.
    pub fn parent(&self, level: u64) -> Self {
        let lsb = lsb_for_level(level);
        CellID((self.0 & (-(lsb as i64)) as u64) | lsb)
    }

    /// immediate_parent is cheaper than Parent, but assumes !ci.isFace().
    pub fn immediate_parent(&self) -> Self {
        let nlsb = self.lsb() << 2;
        CellID((self.0 & (-(nlsb as i64)) as u64) | nlsb)
    }

    //TODO private
    /// is_face returns whether this is a top-level (face) cell.
    pub fn is_face(&self) -> bool {
        (self.0 & (lsb_for_level(0) - 1)) == 0
    }

    //TODO private
    /// lsb returns the least significant bit that is set.
    pub fn lsb(&self) -> u64 {
        self.0 & (-(self.0 as i64) as u64)
    }

    /// children returns the four immediate children of this cell.
    /// If ci is a leaf cell, it returns four identical cells that are not the children.
    pub fn children(&self) -> [CellID; 4] {
        let mut lsb = self.lsb();
        let ch0 = self.0 - lsb + (lsb >> 2);
        lsb >>= 1;
        let ch1 = ch0 + lsb;
        let ch2 = ch1 + lsb;
        let ch3 = ch2 + lsb;

        [CellID(ch0), CellID(ch1), CellID(ch2), CellID(ch3)]
    }

    //TODO private
    /// face_ij_orientation uses the global lookupIJ table to unfiddle the bits of ci.
    pub fn face_ij_orientation(&self) -> (u8, i32, i32, u8) {
        let f = self.face();
        let mut i = 0i32;
        let mut j = 0i32;
        let mut orientation = (f & SWAP_MASK) as u64;
        let mut nbits = MAX_LEVEL - 7 * LOOKUP_BITS;

        let mut k = 7;
        loop {
            orientation += (((self.0 >> (k * 2 * LOOKUP_BITS + 1)) &
                             ((1 << (2 * nbits)) - 1)) as u64) << 2;

            orientation = LOOKUP_IJ[orientation as usize] as u64;
            i += ((orientation as i32) >> (LOOKUP_BITS + 2)) << (k * LOOKUP_BITS);
            j += (((orientation as i32) >> 2) & ((1 << LOOKUP_BITS) - 1)) << (k * LOOKUP_BITS);
            orientation &= (SWAP_MASK | INVERT_MASK) as u64;
            nbits = LOOKUP_BITS; // following iterations

            if k == 0 {
                break;
            }
            k -= 1;
        }

        if self.lsb() & 0x1111111111111110 != 0 {
            orientation ^= SWAP_MASK as u64;
        }
        return (f, i, j, orientation as u8);
    }

    /// edge_neighbors returns the four cells that are adjacent across the cell's four edges.
    /// Edges 0, 1, 2, 3 are in the down, right, up, left directions in the face space.
    /// All neighbors are guaranteed to be distinct.
    pub fn edge_neighbors(&self) -> [CellID; 4] {
        let level = self.level();
        let size = size_ij(level) as i32;
        let (f, i, j, _) = self.face_ij_orientation();

        [CellID::from_face_ij_wrap(f, i, j - size).parent(level),
         CellID::from_face_ij_wrap(f, i + size, j).parent(level),
         CellID::from_face_ij_wrap(f, i, j + size).parent(level),
         CellID::from_face_ij_wrap(f, i - size, j).parent(level)]
    }

    /// vertex_neighbors returns the neighboring cellIDs with vertex closest to this cell at the given level.
    /// (Normally there are four neighbors, but the closest vertex may only have three neighbors if it is one of
    /// the 8 cube vertices.)
    pub fn vertex_neighbors(&self, level: u64) -> Vec<CellID> {
        let half_size = size_ij(level + 1) as i32;
        let size = half_size << 1;
        let (f, i, j, _) = self.face_ij_orientation();

        let (isame, ioffset) = if i & half_size != 0 {
            (i + size < MAX_SIZE_I32, size)
        } else {
            (i - size >= 0, -size)
        };
        let (jsame, joffset) = if j & half_size != 0 {
            (j + size < MAX_SIZE_I32, size)
        } else {
            (j - size >= 0, -size)
        };

        let mut results = Vec::with_capacity(4);
        results.push(self.parent(level));
        results.push(CellID::from_face_ij_same(f, i + ioffset, j, isame).parent(level));
        results.push(CellID::from_face_ij_same(f, i, j + joffset, jsame).parent(level));
        if isame || jsame {
            results.push(CellID::from_face_ij_same(f, i + ioffset, j + joffset, isame && jsame)
                             .parent(level));
        }
        results
    }

    /// all_neighbors returns all neighbors of this cell at the given level. Two
    /// cells X and Y are neighbors if their boundaries intersect but their
    /// interiors do not. In particular, two cells that intersect at a single
    /// point are neighbors. Note that for cells adjacent to a face vertex, the
    /// same neighbor may be returned more than once. There could be up to eight
    /// neighbors including the diagonal ones that share the vertex.
    ///
    /// This requires level >= ci.Level().
    pub fn all_neighbors(&self, level: u64) -> Vec<CellID> {
        let mut neighbors = Vec::new();

        let (face, mut i, mut j, _) = self.face_ij_orientation();

        let size = size_ij(self.level()) as i32;
        i &= -size;
        j &= -size;

        let nbr_size = size_ij(level) as i32;

        let mut k = -nbr_size;
        loop {
            let same_face = if k < 0 {
                (j + k) >= 0
            } else if k >= size {
                (j + k) < MAX_SIZE_I32
            } else {
                neighbors.push(CellID::from_face_ij_same(face, i + k, j - nbr_size, j - size >= 0)
                                   .parent(level));
                neighbors.push(CellID::from_face_ij_same(face,
                                                         i + k,
                                                         j + size,
                                                         j + size < MAX_SIZE_I32)
                                       .parent(level));
                true
            };

            neighbors.push(CellID::from_face_ij_same(face,
                                                     i - nbr_size,
                                                     j + k,
                                                     same_face && i - size >= 0)
                                   .parent(level));
            neighbors.push(CellID::from_face_ij_same(face,
                                                     i + size,
                                                     j + k,
                                                     same_face && i + size < MAX_SIZE_I32)
                                   .parent(level));

            if k >= size {
                break;
            }
            k += nbr_size;
        }

        neighbors
    }

    /// range_min returns the minimum CellID that is contained within this cell.
    pub fn range_min(&self) -> Self {
        CellID(self.0 - (self.lsb() - 1))
    }

    /// range_max returns the maximum CellID that is contained within this cell.
    pub fn range_max(&self) -> Self {
        CellID(self.0 + (self.lsb() - 1))
    }

    /// contains returns true iff the CellID contains oci.
    pub fn contains(&self, other: &CellID) -> bool {
        &self.range_min() <= other && other <= &self.range_max()
    }

    /// intersects returns true iff the CellID intersects oci.
    pub fn intersects(&self, other: &CellID) -> bool {
        other.range_min() <= self.range_max() && other.range_max() >= self.range_min()
    }

    /// face_siti returns the Face/Si/Ti coordinates of the center of the cell.
    fn face_siti(&self) -> (u8, i32, i32) {
        let (face, i, j, _) = self.face_ij_orientation();
        let delta = if self.is_leaf() {
            1
        } else if (i ^ (self.0 as i32 >> 2)) & 1 != 0 {
            2
        } else {
            0
        };
        (face, 2 * i + delta, 2 * j + delta)
    }

    //TODO private
    pub fn raw_point(&self) -> Vector {
        let (face, si, ti) = self.face_siti();
        face_uv_to_xyz(face,
                       st_to_uv(siti_to_st(si as u64)),
                       st_to_uv(siti_to_st(ti as u64)))
    }

    /// child_begin returns the first child in a traversal of the children of this cell, in Hilbert curve order.
    pub fn child_begin(&self) -> Self {
        let ol = self.lsb();
        CellID(self.0 - ol + (ol >> 2))
    }

    /// child_begin_at_level returns the first cell in a traversal of children a given level deeper than this cell, in
    /// Hilbert curve order. The given level must be no smaller than the cell's level.
    /// See ChildBegin for example use.
    pub fn child_begin_at_level(&self, level: u64) -> Self {
        assert!(self.level() <= level);

        CellID(self.0 - self.lsb() + lsb_for_level(level))
    }

    /// child_end returns the first cell after a traversal of the children of this cell in Hilbert curve order.
    /// The returned cell may be invalid.
    pub fn child_end(&self) -> Self {
        let ol = self.lsb();
        CellID(self.0 + ol + (ol >> 2))
    }

    /// child_end_at_level returns the first cell after the last child in a traversal of children a given level deeper
    /// than this cell, in Hilbert curve order.
    /// The given level must be no smaller than the cell's level.
    /// The returned cell may be invalid.
    pub fn child_end_at_level(&self, level: u64) -> Self {
        assert!(self.level() <= level);

        CellID(self.0 + self.lsb() + lsb_for_level(level))
    }

    /// next returns the next cell along the Hilbert curve.
    /// This is expected to be used with child_begin and child_end,
    /// or child_begin_at_level and child_end_at_level.
    pub fn next(&self) -> Self {
        CellID(self.0.wrapping_add(self.lsb() << 1))
    }

    /// prev returns the previous cell along the Hilbert curve.
    pub fn prev(&self) -> Self {
        CellID(self.0.wrapping_sub(self.lsb() << 1))
    }

    /// next_wrap returns the next cell along the Hilbert curve, wrapping from last to
    /// first as necessary. This should not be used with ChildBegin and ChildEnd.
    pub fn next_wrap(&self) -> Self {
        let next = self.next();
        if next.0 < WRAP_OFFSET {
            next
        } else {
            CellID(next.0.wrapping_sub(WRAP_OFFSET))
        }
    }

    /// prev_wrap returns the previous cell along the Hilbert curve, wrapping around from
    /// first to last as necessary. This should not be used with ChildBegin and ChildEnd.
    pub fn prev_wrap(&self) -> Self {
        let prev = self.prev();
        if prev.0 < WRAP_OFFSET {
            prev
        } else {
            CellID(prev.0.wrapping_add(WRAP_OFFSET))
        }
    }

    /// advance_wrap advances or retreats the indicated number of steps along the
    /// Hilbert curve at the current level and returns the new position. The
    /// position wraps between the first and last faces as necessary.
    pub fn advance_wrap(&self, mut steps: i64) -> Self {
        if steps == 0 {
            return self.clone();
        }

        let shift = 2 * (MAX_LEVEL - self.level()) + 1;
        if steps < 0 {
            let min = -((self.0 >> shift) as i64);
            if steps < min {
                let wrap = (WRAP_OFFSET >> shift) as i64;
                steps %= wrap;
                if steps < min {
                    steps += wrap;
                }
            }
        } else {
            let max = ((WRAP_OFFSET - self.0) >> shift) as i64;
            if steps > max {
                let wrap = (WRAP_OFFSET >> shift) as i64;
                steps %= wrap;
                if steps > max {
                    steps -= wrap;
                }
            }
        }
        CellID(self.0.wrapping_add((steps as u64) << shift))
    }

    // TODO: the methods below are not exported yet.  Settle on the entire API design
    // before doing this.  Do we want to mirror the C++ one as closely as possible?

    //TODO private
    /// distance_from_begin returns the number of steps that this cell is from the first
    /// node in the S2 heirarchy at our level. (i.e., FromFace(0).ChildBeginAtLevel(ci.Level())).
    /// The return value is always non-negative.
    pub fn distance_from_begin(&self) -> i64 {
        (self.0 >> (2 * (MAX_LEVEL - self.level()) + 1)) as i64
    }

    /// common_ancestor_level returns the level of the common ancestor of the two S2 CellIDs.
    pub fn common_ancestor_level(&self, other: &Self) -> Option<u64> {
        let mut bits = self.0 ^ other.0;

        if bits < self.lsb() {
            bits = self.lsb();
        }
        if bits < other.lsb() {
            bits = other.lsb();
        }

        let msb_pos = find_msb_set_nonzero64(bits) as u64;
        if msb_pos > 60 {
            None
        } else {
            Some((60 - msb_pos) >> 1)
        }
    }

    /// advance advances or retreats the indicated number of steps along the
    /// Hilbert curve at the current level, and returns the new position. The
    /// position is never advanced past End() or before Begin().
    pub fn advance(&self, mut steps: i64) -> Self {
        if steps == 0 {
            return self.clone();
        }

        let step_shift = (2 * (MAX_LEVEL - self.level()) + 1) as u64;
        if steps < 0 {
            let min_steps = -((self.0 >> step_shift) as i64);
            if steps < min_steps {
                steps = min_steps;
            }
        } else {
            let max_steps = ((WRAP_OFFSET + self.lsb() - self.0) >> step_shift) as i64;
            if steps > max_steps {
                steps = max_steps;
            }
        }
        CellID(self.0.wrapping_add((steps << step_shift) as u64))
    }

    /// center_st return the center of the CellID in (s,t)-space.
    fn center_st(&self) -> r2::Point {
        let (_, si, ti) = self.face_siti();
        r2::Point {
            x: siti_to_st(si as u64),
            y: siti_to_st(ti as u64),
        }
    }

    /// size_st returns the edge length of this CellID in (s,t)-space at the given level.
    fn size_st(&self, level: u64) -> f64 {
        ij_to_stmin(size_ij(level) as i32)
    }

    //TODO private
    /// bound_st returns the bound of this CellID in (s,t)-space.
    pub fn bound_st(&self) -> r2::Rect {
        let s = self.size_st(self.level());
        r2::Rect::from_center_size(&self.center_st(), &r2::Point { x: s, y: s })
    }

    //TODO private
    /// center_uv returns the center of this CellID in (u,v)-space. Note that
    /// the center of the cell is defined as the point at which it is recursively
    /// subdivided into four children; in general, it is not at the midpoint of
    /// the (u,v) rectangle covered by the cell.
    pub fn center_uv(&self) -> r2::Point {
        let (_, si, ti) = self.face_siti();
        r2::Point {
            x: st_to_uv(siti_to_st(si as u64)),
            y: st_to_uv(siti_to_st(ti as u64)),
        }
    }

    /// bound_uv returns the bound of this CellID in (u,v)-space.
    pub fn bound_uv(&self) -> r2::Rect {
        let (_, i, j, _) = self.face_ij_orientation();
        ij_level_to_bound_uv(i, j, self.level())
    }

    //TODO private
    pub fn expand_by_distance_uv(uv: r2::Rect, distance: Angle) -> r2::Rect {
        let max_u = uv.x.lo.abs().max(uv.x.hi.abs());
        let max_v = uv.y.lo.abs().max(uv.y.hi.abs());

        let sin_dist = distance.0.sin();
        r2::Rect {
            x: r1::Interval {
                lo: expand_endpoint(uv.x.lo, max_v, -sin_dist),
                hi: expand_endpoint(uv.x.hi, max_v, sin_dist),
            },
            y: r1::Interval {
                lo: expand_endpoint(uv.y.lo, max_u, -sin_dist),
                hi: expand_endpoint(uv.y.hi, max_u, sin_dist),
            },
        }
    }

    /// max_tile returns the largest cell with the same range_min such that
    /// range_max < limit.range_min. It returns limit if no such cell exists.
    /// This method can be used to generate a small set of CellIDs that covers
    /// a given range (a tiling). This example shows how to generate a tiling
    /// for a semi-open range of leaf cells [start, limit):
    ///
    /// for id := start.MaxTile(limit); id != limit; id = id.Next().MaxTile(limit)) { ... }
    ///
    /// Note that in general the cells in the tiling will be of different sizes;
    /// they gradually get larger (near the middle of the range) and then
    /// gradually get smaller as limit is approached.
    pub fn max_tile(&self, limit: &Self) -> Self {
        let mut s = self.clone();
        let start = s.range_min();
        if start >= limit.range_min() {
            return limit.clone();
        }
        if s.range_max() > *limit {
            loop {
                s = s.children()[0].clone();
                if s.range_max() < *limit {
                    break;
                }
            }
            return s.clone();
        }

        while !s.is_face() {
            let parent = self.immediate_parent();
            if parent.range_min() != start || parent.range_max() >= *limit {
                break;
            }
            s = parent;
        }
        s.clone()
    }

    //TODO private
    /// center_face_siti returns the (face, si, ti) coordinates of the center of the cell.
    /// Note that although (si,ti) coordinates span the range [0,2**31] in general,
    /// the cell center coordinates are always in the range [1,2**31-1] and
    /// therefore can be represented using a signed 32-bit integer.
    pub fn center_face_siti(&self) -> (u8, i32, i32) {
        let (face, i, j, _) = self.face_ij_orientation();
        let delta = if self.is_leaf() {
            1
        } else if ((i as i64) ^ (self.0 as i64 >> 2)) & 1 == 1 {
            2
        } else {
            0
        };

        (face, 2 * i + delta, 2 * j + delta)
    }
}

/// expand_endpoint returns a new u-coordinate u' such that the distance from the
/// line u=u' to the given edge (u,v0)-(u,v1) is exactly the given distance
/// (which is specified as the sine of the angle corresponding to the distance).
fn expand_endpoint(u: f64, max_v: f64, sin_dist: f64) -> f64 {
    let sin_u_shift = sin_dist * ((1. + u * u + max_v * max_v) / (1. + u * u)).sqrt();
    let cos_u_shift = (1. - sin_u_shift * sin_u_shift).sqrt();
    (cos_u_shift * u + sin_u_shift) / (cos_u_shift - sin_u_shift * u)
}

/// ij_tostmin converts the i- or j-index of a leaf cell to the minimum corresponding
/// s- or t-value contained by that cell. The argument must be in the range
/// [0..2**30], i.e. up to one position beyond the normal range of valid leaf
/// cell indices.
fn ij_to_stmin(i: i32) -> f64 {
    return (i as f64) / (MAX_SIZE as f64);
}

/// st_to_ij converts value in ST coordinates to a value in IJ coordinates.
fn st_to_ij(s: f64) -> i32 {
    clamp(((MAX_SIZE as f64 * s).floor() as i32), 0, MAX_SIZE_I32 - 1)
}

impl std::fmt::Debug for CellID {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}/", self.face())?;
        for level in 1..(self.level() + 1) {
            write!(f, "{}", self.child_position(level).to_string())?;
        }
        Ok(())
    }
}

impl From<CellID> for Point {
    /// returns the center of the s2 cell on the sphere as a Point.
    /// The maximum directional error in Point (compared to the exact
    /// mathematical result) is 1.5 * dblEpsilon radians, and the maximum length
    /// error is 2 * dblEpsilon (the same as Normalize).
    fn from(id: CellID) -> Self {
        Point::from(&id)
    }
}
impl<'a> From<&'a CellID> for Point {
    /// cellIDFromPoint returns a leaf cell containing point p. Usually there is
    /// exactly one such cell, but for points along the edge of a cell, any
    /// adjacent cell may be (deterministically) chosen. This is because
    /// s2.CellIDs are considered to be closed sets. The returned cell will
    /// always contain the given point, i.e.
    ///
    /// Cell::from(&p).contains_point(&p)
    ///
    /// is always true.
    fn from(id: &'a CellID) -> Self {
        Point(id.raw_point().normalize())
    }
}

impl From<CellID> for LatLng {
    /// returns the center of the s2 cell on the sphere as a LatLng.
    fn from(id: CellID) -> Self {
        LatLng::from(&id)
    }
}
impl<'a> From<&'a CellID> for LatLng {
    fn from(id: &'a CellID) -> Self {
        LatLng::from(Point::from(id))
    }
}


impl From<LatLng> for CellID {
    fn from(ll: LatLng) -> Self {
        let p: Point = ll.into();
        Self::from(p)
    }
}

impl<'a> From<&'a Point> for CellID {
    fn from(p: &'a Point) -> Self {
        let (f, u, v) = xyz_to_face_uv(&p.0);
        let i = st_to_ij(uv_to_st(u));
        let j = st_to_ij(uv_to_st(v));
        return CellID::from_face_ij(f, i, j);
    }
}
impl From<Point> for CellID {
    fn from(p: Point) -> Self {
        CellID::from(&p)
    }
}

pub struct CellIDIter {
    cur: CellID,
    end: CellID,
}

impl Iterator for CellIDIter {
    type Item = CellID;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur == self.end {
            None
        } else {
            let res = Some(self.cur.clone());
            self.cur = self.cur.next();
            res
        }
    }
}

impl CellID {
    pub fn child_iter(&self) -> CellIDIter {
        CellIDIter {
            cur: self.child_begin(),
            end: self.child_end(),
        }
    }

    pub fn child_iter_at_level(&self, level: u64) -> CellIDIter {
        CellIDIter {
            cur: self.child_begin_at_level(level),
            end: self.child_end_at_level(level),
        }
    }
}

//TODO private
pub const IJ_TO_POS: [[u8; 4]; 4] = [[0, 1, 3, 2], [0, 3, 1, 2], [2, 3, 1, 0], [2, 1, 3, 0]];
pub const POS_TO_IJ: [[u8; 4]; 4] = [[0, 1, 3, 2], [0, 2, 3, 1], [3, 2, 0, 1], [3, 1, 0, 2]];
pub const POS_TO_ORIENTATION: [u8; 4] = [SWAP_MASK, 0, 0, INVERT_MASK | SWAP_MASK];

lazy_static! {
    static ref LOOKUP_TBL: [Vec<u64>; 2] = {
        let size = 1 << (2*LOOKUP_BITS + 2);
        let mut lookup_pos = Vec::new();
        let mut lookup_ij = Vec::new();
        lookup_pos.resize(size, 0);
        lookup_ij.resize(size, 0);

        init_lookup_cell(0, 0, 0, 0, 0, 0, &mut lookup_pos, &mut lookup_ij);
        init_lookup_cell(0, 0, 0, SWAP_MASK, 0, SWAP_MASK, &mut lookup_pos, &mut lookup_ij);
        init_lookup_cell(0, 0, 0, INVERT_MASK, 0, INVERT_MASK, &mut lookup_pos, &mut lookup_ij);
        init_lookup_cell(0, 0, 0, SWAP_MASK|INVERT_MASK, 0, SWAP_MASK|INVERT_MASK, &mut lookup_pos, &mut lookup_ij);
        [lookup_pos, lookup_ij]
    };

    static ref LOOKUP_POS: &'static [u64] = {
        LOOKUP_TBL[0].as_slice()
    };

    static ref LOOKUP_IJ: &'static [u64] = {
        LOOKUP_TBL[1].as_slice()
    };
}

/// init_lookup_cell initializes the lookupIJ table at init time.
fn init_lookup_cell(level: u64,
                    i: i32,
                    j: i32,
                    orig_orientation: u8,
                    pos: usize,
                    orientation: u8,
                    lookup_pos: &mut [u64],
                    lookup_ij: &mut [u64]) {
    if level == LOOKUP_BITS {
        let ij = ((i << LOOKUP_BITS) + j) as usize;
        lookup_pos[(ij << 2) + orig_orientation as usize] = (pos << 2) as u64 + orientation as u64;
        lookup_ij[(pos << 2) + orig_orientation as usize] = (ij << 2) as u64 + orientation as u64;
        return;
    }

    let r = &POS_TO_IJ[orientation as usize];
    for idx in 0..4 {
        init_lookup_cell(level + 1,
                         (i << 1) + (r[idx] >> 1) as i32,
                         (j << 1) + (r[idx] & 1) as i32,
                         orig_orientation,
                         (pos << 2) + idx,
                         orientation ^ POS_TO_ORIENTATION[idx],
                         lookup_pos,
                         lookup_ij)
    }
}

/// ij_level_to_bound_uv returns the bounds in (u,v)-space for the cell at the given
/// level containing the leaf cell with the given (i,j)-coordinates.
pub fn ij_level_to_bound_uv(i: i32, j: i32, level: u64) -> r2::Rect {
    let cell_size = size_ij(level) as i32;
    let x_lo = i & -cell_size;
    let y_lo = j & -cell_size;

    r2::Rect {
        x: r1::Interval {
            lo: st_to_uv(ij_to_stmin(x_lo)),
            hi: st_to_uv(ij_to_stmin(x_lo + cell_size)),
        },
        y: r1::Interval {
            lo: st_to_uv(ij_to_stmin(y_lo)),
            hi: st_to_uv(ij_to_stmin(y_lo + cell_size)),
        },
    }
}

fn find_msb_set_nonzero64(bits: u64) -> u32 {
    63 - bits.leading_zeros()
}
//TODO private
pub fn find_lsb_set_nonzero64(bits: u64) -> u32 {
    bits.trailing_zeros()
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use rand::Rng;
    use s1;
    use s2::random;

    #[test]
    fn test_cellid_from_face() {
        for face in 0..6 {
            let fpl = CellID::from_face_pos_level(face, 0, 0);
            let f = CellID::from_face(face);
            assert_eq!(fpl, f);
        }
    }

    #[test]
    fn test_cellid_parent_child_relationships() {
        let ci = CellID::from_face_pos_level(3, 0x12345678, MAX_LEVEL - 4);
        assert!(ci.is_valid());
        assert_eq!(ci.face(), 3);
        assert_eq!(ci.pos(), 0x12345700);
        assert_eq!(ci.level(), 26);
        assert_eq!(ci.is_leaf(), false);

        assert_eq!(ci.child_begin_at_level(ci.level() + 2).pos(), 0x12345610);
        assert_eq!(ci.child_begin().pos(), 0x12345640);
        assert_eq!(ci.children()[0].pos(), 0x12345640);
        assert_eq!(ci.immediate_parent().pos(), 0x12345400);
        assert_eq!(ci.parent(ci.level() - 2).pos(), 0x12345000);

        assert!(ci.child_begin() < ci);
        assert!(ci.child_end() > ci);
        assert_eq!(ci.child_end(), ci.child_begin().next().next().next().next());

        assert_eq!(ci.range_min(), ci.child_begin_at_level(MAX_LEVEL));
        assert_eq!(ci.range_max().next(), ci.child_end_at_level(MAX_LEVEL));
    }

    fn test_cellid_containment_case(x: &CellID,
                                    y: &CellID,
                                    x_contains_y: bool,
                                    y_contains_x: bool,
                                    x_intersects_y: bool) {
        assert_eq!(x.contains(y), x_contains_y);
        assert_eq!(y.contains(x), y_contains_x);
        assert_eq!(x.intersects(y), x_intersects_y);
        assert_eq!(y.intersects(x), x_intersects_y);
    }

    #[test]
    fn test_cellid_containment() {
        let a = CellID(0x80855c0000000000); // Pittsburg
        let b = CellID(0x80855d0000000000); // child of a
        let c = CellID(0x80855dc000000000); // child of b
        let d = CellID(0x8085630000000000); // part of Pittsburg disjoint from a

        test_cellid_containment_case(&a, &a, true, true, true);
        test_cellid_containment_case(&a, &b, true, false, true);
        test_cellid_containment_case(&a, &c, true, false, true);
        test_cellid_containment_case(&a, &d, false, false, false);
        test_cellid_containment_case(&b, &b, true, true, true);
        test_cellid_containment_case(&b, &c, true, false, true);
        test_cellid_containment_case(&b, &d, false, false, false);
        test_cellid_containment_case(&c, &c, true, true, true);
        test_cellid_containment_case(&c, &d, false, false, false);
        test_cellid_containment_case(&d, &d, true, true, true);

        // TODO(dsymonds): Test Contains, Intersects better, such as with adjacent cells.
    }

    #[test]
    fn test_cellid_string() {
        let ci = CellID(0xbb04000000000000);
        assert_eq!(format!("{:?}", ci), "5/31200");
    }

    fn test_cellid_latlng_case(ci: CellID, lat: f64, lng: f64) {
        let ll = LatLng {
            lat: s1::angle::Deg(lat).into(),
            lng: s1::angle::Deg(lng).into(),
        };
        //TODO
        let l2: LatLng = ci.clone().into();

        let distance = ll.distance(&l2);
        assert!(distance < s1::angle::Deg(1.0e-9).into());

        let ci2: CellID = ll.into();
        assert_eq!(ci, ci2);
    }

    #[test]
    fn test_cellid_latlng() {
        test_cellid_latlng_case(CellID(0x47a1cbd595522b39), 49.703498679, 11.770681595);
        test_cellid_latlng_case(CellID(0x46525318b63be0f9), 55.685376759, 12.588490937);
        test_cellid_latlng_case(CellID(0x52b30b71698e729d), 45.486546517, -93.449700022);
        test_cellid_latlng_case(CellID(0x46ed8886cfadda85), 58.299984854, 23.049300056);
        test_cellid_latlng_case(CellID(0x3663f18a24cbe857), 34.364439040, 108.330699969);
        test_cellid_latlng_case(CellID(0x10a06c0a948cf5d), -30.694551352, -30.048758753);
        test_cellid_latlng_case(CellID(0x2b2bfd076787c5df), -25.285264027, 133.823116966);
        test_cellid_latlng_case(CellID(0xb09dff882a7809e1), -75.000000031, 0.000000133);
        test_cellid_latlng_case(CellID(0x94daa3d000000001), -24.694439215, -47.537363213);
        test_cellid_latlng_case(CellID(0x87a1000000000001), 38.899730392, -99.901813021);
        test_cellid_latlng_case(CellID(0x4fc76d5000000001), 81.647200334, -55.631712940);
        test_cellid_latlng_case(CellID(0x3b00955555555555), 10.050986518, 78.293170610);
        test_cellid_latlng_case(CellID(0x1dcc469991555555), -34.055420593, 18.551140038);
        test_cellid_latlng_case(CellID(0xb112966aaaaaaaab), -69.219262171, 49.670072392);
    }

    #[test]
    fn test_cellid_edge_neighbors() {
        let faces = [5, 3, 2, 0];

        for (i, nbr) in CellID::from_face_ij(1, 0, 0)
                .parent(0)
                .edge_neighbors()
                .iter()
                .enumerate() {
            assert!(nbr.is_face());
            assert_eq!(nbr.face(), faces[i]);
        }

        let max_ij = MAX_SIZE_I32 - 1;
        for level in 1..(MAX_LEVEL + 1) {
            let id = CellID::from_face_ij(1, 0, 0).parent(level);
            let level_size_ij = size_ij(level) as i32;
            let want = [CellID::from_face_ij(5, max_ij, max_ij).parent(level),
                        CellID::from_face_ij(1, level_size_ij, 0).parent(level),
                        CellID::from_face_ij(1, 0, level_size_ij).parent(level),
                        CellID::from_face_ij(0, max_ij, 0).parent(level)];

            assert_eq!(want, id.edge_neighbors());
        }
    }

    #[test]
    fn test_cellid_vertex_neighbors() {
        let id: CellID = Point(Vector {
                                   x: 0.,
                                   y: 0.,
                                   z: 1.,
                               })
                .into();

        let mut neighbors = id.vertex_neighbors(5);
        neighbors.sort();

        for (n, nbr) in neighbors.iter().enumerate() {
            let mut i = 1 << 29;
            let mut j = 1 << 29;
            if n < 2 {
                i -= 1;
            }
            if n == 0 || n == 3 {
                j -= 1;
            }
            assert_eq!(*nbr, CellID::from_face_ij(2, i, j).parent(5));
        }

        let id2 = CellID::from_face_pos_level(0, 0, MAX_LEVEL);
        let mut neighbors2 = id2.vertex_neighbors(0);
        neighbors2.sort();
        assert!(neighbors2.len() == 3);
        assert_eq!(neighbors2[0], CellID::from_face(0));
        assert_eq!(neighbors2[1], CellID::from_face(4));
    }

    #[test]
    fn test_cellid_all_neighbors() {
        let mut rng = random::rng();
        //XXX TODO FIX to 1000
        for _ in 0..10 {
            let mut id = random::cellid(&mut rng);
            if id.is_leaf() {
                id = id.immediate_parent();
            }

            let mut max_diff = MAX_LEVEL - id.level() - 1;
            if max_diff > 6 {
                max_diff = 6;
            }
            let level = match max_diff {
                0 => id.level(),
                _ => id.level() + rng.gen_range(0, max_diff),
            };

            let mut want = Vec::new();
            let mut all = id.all_neighbors(level);

            let end = id.child_end_at_level(level + 1);
            let mut c = id.child_begin_at_level(level + 1);
            while c != end {
                all.push(c.immediate_parent());
                want.extend_from_slice(&c.vertex_neighbors(level));
                c = c.next();
            }

            all.sort();
            all.dedup();
            want.sort();
            want.dedup();

            assert_eq!(all, want);
        }
    }

    fn test_cellid_tokens_case(s: &str, id: CellID) {
        assert_eq!(CellID::from_token(s), id);
        assert_eq!(s, id.to_token());
    }

    #[test]
    fn test_cellid_tokens_nominal() {
        test_cellid_tokens_case("1", CellID(0x1000000000000000));
        test_cellid_tokens_case("3", CellID(0x3000000000000000));
        test_cellid_tokens_case("14", CellID(0x1400000000000000));
        test_cellid_tokens_case("41", CellID(0x4100000000000000));
        test_cellid_tokens_case("094", CellID(0x0940000000000000));
        test_cellid_tokens_case("537", CellID(0x5370000000000000));
        test_cellid_tokens_case("3fec", CellID(0x3fec000000000000));
        test_cellid_tokens_case("72f3", CellID(0x72f3000000000000));
        test_cellid_tokens_case("52b8c", CellID(0x52b8c00000000000));
        test_cellid_tokens_case("990ed", CellID(0x990ed00000000000));
        test_cellid_tokens_case("4476dc", CellID(0x4476dc0000000000));
        test_cellid_tokens_case("2a724f", CellID(0x2a724f0000000000));
        test_cellid_tokens_case("7d4afc4", CellID(0x7d4afc4000000000));
        test_cellid_tokens_case("b675785", CellID(0xb675785000000000));
        test_cellid_tokens_case("40cd6124", CellID(0x40cd612400000000));
        test_cellid_tokens_case("3ba32f81", CellID(0x3ba32f8100000000));
        test_cellid_tokens_case("08f569b5c", CellID(0x08f569b5c0000000));
        test_cellid_tokens_case("385327157", CellID(0x3853271570000000));
        test_cellid_tokens_case("166c4d1954", CellID(0x166c4d1954000000));
        test_cellid_tokens_case("96f48d8c39", CellID(0x96f48d8c39000000));
        test_cellid_tokens_case("0bca3c7f74c", CellID(0x0bca3c7f74c00000));
        test_cellid_tokens_case("1ae3619d12f", CellID(0x1ae3619d12f00000));
        test_cellid_tokens_case("07a77802a3fc", CellID(0x07a77802a3fc0000));
        test_cellid_tokens_case("4e7887ec1801", CellID(0x4e7887ec18010000));
        test_cellid_tokens_case("4adad7ae74124", CellID(0x4adad7ae74124000));
        test_cellid_tokens_case("90aba04afe0c5", CellID(0x90aba04afe0c5000));
        test_cellid_tokens_case("8ffc3f02af305c", CellID(0x8ffc3f02af305c00));
        test_cellid_tokens_case("6fa47550938183", CellID(0x6fa4755093818300));
        test_cellid_tokens_case("aa80a565df5e7fc", CellID(0xaa80a565df5e7fc0));
        test_cellid_tokens_case("01614b5e968e121", CellID(0x01614b5e968e1210));
        test_cellid_tokens_case("aa05238e7bd3ee7c", CellID(0xaa05238e7bd3ee7c));
        test_cellid_tokens_case("48a23db9c2963e5b", CellID(0x48a23db9c2963e5b));
    }

    #[test]
    fn test_cellid_tokens_error_case() {
        assert_eq!("X", CellID(0).to_token());
        assert_eq!(CellID(0), CellID::from_token("X"));

        assert_eq!(CellID(0), CellID::from_token("876b e99"));
        assert_eq!(CellID(0), CellID::from_token("876bee99\n"));
        assert_eq!(CellID(0), CellID::from_token("876[ee99"));
        assert_eq!(CellID(0), CellID::from_token(" 876bee99"));
    }

    //TODO: duplicated
    fn test_approx_eq(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-14);
    }

    macro_rules! P {
        ($x: expr, $y: expr) => {
            r2::Point::new($x as f64, $y as f64)
        }
    }

    fn test_ij_level_to_bound_uv_case(i: i32, j: i32, level: u64, points: &[r2::Point]) {
        let uv = ij_level_to_bound_uv(i, j, level);
        let want = r2::Rect::from_points(points);

        test_approx_eq(uv.x.lo, want.x.lo);
        test_approx_eq(uv.x.hi, want.x.hi);

        test_approx_eq(uv.y.lo, want.y.lo);
        test_approx_eq(uv.y.hi, want.y.hi);
    }

    const MAX_IJ: i32 = MAX_SIZE_I32 - 1;

    #[test]
    fn test_ij_level_to_bound_uv() {
        test_ij_level_to_bound_uv_case(-1, -1, 0, &[P!(-5., -5.), P!(-1., -1.)]);
        test_ij_level_to_bound_uv_case(-1 * MAX_IJ, -1 * MAX_IJ, 0, &[P!(-5., -5.), P!(-1., -1.)]);
        test_ij_level_to_bound_uv_case(-1,
                                       -1,
                                       MAX_LEVEL,
                                       &[P!(-1.0000000024835267, -1.0000000024835267),
                                         P!(-1., -1.)]);

        //XXX
        // test_ij_level_to_bound_uv_case(0, 0, MAX_LEVEL + 1, &[P!(-1., -1.), P!(-1., -1.)]);

        // Minimum i,j at different levels
        test_ij_level_to_bound_uv_case(0, 0, 0, &[P!(-1., -1.), P!(1., 1.)]);
        test_ij_level_to_bound_uv_case(0,
                                       0,
                                       MAX_LEVEL / 2,
                                       &[P!(-1, -1),
                                         P!(-0.999918621033430099, -0.999918621033430099)]);
        test_ij_level_to_bound_uv_case(0,
                                       0,
                                       MAX_LEVEL,
                                       &[P!(-1, -1),
                                         P!(-0.999999997516473060, -0.999999997516473060)]);
        test_ij_level_to_bound_uv_case(1, 1, 0, &[P!(-1., -1.), P!(1., 1.)]);
        test_ij_level_to_bound_uv_case(1,
                                       1,
                                       MAX_LEVEL / 2,
                                       &[P!(-1, -1),
                                         P!(-0.999918621033430099, -0.999918621033430099)]);

        // Just a hair off the outer bounds at different levels.
        test_ij_level_to_bound_uv_case(1,
                                       1,
                                       MAX_LEVEL,
                                       &[P!(-0.9999999975164731, -0.9999999975164731),
                                         P!(-0.9999999950329462, -0.9999999950329462)]);

        test_ij_level_to_bound_uv_case(MAX_IJ / 2, MAX_IJ / 2, 0, &[P!(-1, -1), P!(1, 1)]);

        test_ij_level_to_bound_uv_case(MAX_IJ / 2,
                                       MAX_IJ / 2,
                                       MAX_LEVEL / 2,
                                       &[P!(-0.000040691345930099, -0.000040691345930099),
                                         P!(0., 0.)]);

        test_ij_level_to_bound_uv_case(MAX_IJ / 2,
                                       MAX_IJ / 2,
                                       MAX_LEVEL,
                                       &[P!(-0.000000001241763433, -0.000000001241763433),
                                         P!(0., 0.)]);

        // Maximum i, j at different levels.
        test_ij_level_to_bound_uv_case(MAX_IJ, MAX_IJ, 0, &[P!(-1., -1.), P!(1., 1.)]);

        // Center point of the i,j space at different levels.
        test_ij_level_to_bound_uv_case(MAX_IJ,
                                       MAX_IJ,
                                       MAX_LEVEL / 2,
                                       &[P!(0.999918621033430099, 0.999918621033430099),
                                         P!(1., 1.)]);

        test_ij_level_to_bound_uv_case(MAX_IJ,
                                       MAX_IJ,
                                       MAX_LEVEL,
                                       &[P!(0.999999997516473060, 0.999999997516473060),
                                         P!(1., 1.)]);
    }

    fn test_common_ancestor_case(expected: Option<u64>, c1: CellID, c2: CellID) {
        assert_eq!(expected, c1.common_ancestor_level(&c2));
        assert_eq!(expected, c2.common_ancestor_level(&c1));
    }

    #[test]
    fn test_cellid_common_ancestor_level() {
        test_common_ancestor_case(Some(0), CellID::from_face(0), CellID::from_face(0));
        test_common_ancestor_case(Some(30),
                                  CellID::from_face(0).child_begin_at_level(30),
                                  CellID::from_face(0).child_begin_at_level(30));

        // One cell is a descendant of the other.
        test_common_ancestor_case(Some(0),
                                  CellID::from_face(0),
                                  CellID::from_face(0).child_begin_at_level(30));

        test_common_ancestor_case(Some(0),
                                  CellID::from_face(5),
                                  CellID::from_face(5).child_end_at_level(30).prev());

        // No common ancestors.
        test_common_ancestor_case(None, CellID::from_face(0), CellID::from_face(5));

        test_common_ancestor_case(None,
                                  CellID::from_face(2).child_begin_at_level(30),
                                  CellID::from_face(3).child_begin_at_level(20));

        // Common ancestor distinct from both.
        test_common_ancestor_case(Some(8),
                                  CellID::from_face(5)
                                      .child_begin_at_level(9)
                                      .next()
                                      .child_begin_at_level(15),
                                  CellID::from_face(5)
                                      .child_begin_at_level(9)
                                      .child_begin_at_level(20));

        test_common_ancestor_case(Some(1),
                                  CellID::from_face(0)
                                      .child_begin_at_level(2)
                                      .child_begin_at_level(30),
                                  CellID::from_face(0)
                                      .child_begin_at_level(2)
                                      .next()
                                      .child_begin_at_level(5));
    }

    #[test]
    fn test_cellid_distance_to_begin() {
        assert_eq!(6,
                   CellID::from_face(5)
                       .child_end_at_level(0)
                       .distance_from_begin());
        assert_eq!(6 * (1 << (2 * MAX_LEVEL)),
                   CellID::from_face(5)
                       .child_end_at_level(MAX_LEVEL)
                       .distance_from_begin());
        assert_eq!(0,
                   CellID::from_face(0)
                       .child_begin_at_level(0)
                       .distance_from_begin());
        assert_eq!(0,
                   CellID::from_face(0)
                       .child_begin_at_level(MAX_LEVEL)
                       .distance_from_begin());

        let id = CellID::from_face_pos_level(3, 0x12345678, MAX_LEVEL - 4);
        assert_eq!(id,
                   CellID::from_face(0)
                       .child_begin_at_level(id.level())
                       .advance(id.distance_from_begin()));
    }

    #[test]
    fn test_find_msb_set_nonzero64() {
        let mut test_one = 0x8000000000000000;
        let mut test_all = 0xFFFFFFFFFFFFFFFF;
        let mut test_some = 0xFEDCBA9876543210;

        let mut i = 63;
        loop {
            assert_eq!(i, find_msb_set_nonzero64(test_one));
            assert_eq!(i, find_msb_set_nonzero64(test_all));
            assert_eq!(i, find_msb_set_nonzero64(test_some));
            if i == 0 {
                break;
            }
            i -= 1;

            test_one >>= 1;
            test_all >>= 1;
            test_some >>= 1;
        }

        assert_eq!(0, find_msb_set_nonzero64(1));
        //XXX: overflow panics on rust
        // assert_eq!(0, find_msb_set_nonzero64(0));
    }

    #[test]
    fn test_find_lsb_set_nonzero64() {
        let mut test_one = 0x0000000000000001;
        let mut test_all = 0xFFFFFFFFFFFFFFFF;
        let mut test_some = 0x0123456789ABCDEF;

        for i in 0..64 {
            assert_eq!(i, find_lsb_set_nonzero64(test_one));
            assert_eq!(i, find_lsb_set_nonzero64(test_all));
            assert_eq!(i, find_lsb_set_nonzero64(test_some));

            test_one <<= 1;
            test_all <<= 1;
            test_some <<= 1;
        }

        //XXX: nonzero64 is not supposed to work on
        // assert_eq!(0, find_lsb_set_nonzero64(0));
    }

    #[test]
    fn test_cellid_wrapping() {
        let id = CellID::from_face_pos_level(3, 0x12345678, MAX_LEVEL - 4);

        // "test wrap from beginning to end of Hilbert curve",
        assert_eq!(CellID::from_face(5).child_end_at_level(0).prev(),
                   CellID::from_face(0).child_begin_at_level(0).prev_wrap());

        // "smallest end leaf wraps to smallest first leaf using prev_wrap",
        assert_eq!(CellID::from_face_pos_level(5, (!0u64) >> FACE_BITS, MAX_LEVEL),
                   CellID::from_face(0)
                       .child_begin_at_level(MAX_LEVEL)
                       .prev_wrap());
        // "smallest end leaf wraps to smallest first leaf using advance_wrap",
        assert_eq!(CellID::from_face_pos_level(5, (!0u64) >> FACE_BITS, MAX_LEVEL),
                   CellID::from_face(0)
                       .child_begin_at_level(MAX_LEVEL)
                       .advance_wrap(-1));
        // "prev_wrap is the same as advance_wrap(-1)",
        assert_eq!(CellID::from_face(0)
                       .child_begin_at_level(MAX_LEVEL)
                       .advance_wrap(-1),
                   CellID::from_face(0)
                       .child_begin_at_level(MAX_LEVEL)
                       .prev_wrap());
        // "prev + next_wrap stays the same at given level",
        assert_eq!(CellID::from_face(0).child_begin_at_level(4),
                   CellID::from_face(5)
                       .child_end_at_level(4)
                       .prev()
                       .next_wrap());
        // "advance_wrap forward and back stays the same at given level",
        assert_eq!(CellID::from_face(0).child_begin_at_level(4),
                   CellID::from_face(5)
                       .child_end_at_level(4)
                       .advance(-1)
                       .advance_wrap(1));
        // "prev().next_wrap() stays same for first cell at level",
        assert_eq!(CellID::from_face_pos_level(0, 0, MAX_LEVEL),
                   CellID::from_face(5)
                       .child_end_at_level(MAX_LEVEL)
                       .prev()
                       .next_wrap());
        // "advance_wrap forward and back stays same for first cell at level",
        assert_eq!(CellID::from_face_pos_level(0, 0, MAX_LEVEL),
                   CellID::from_face(5)
                       .child_end_at_level(MAX_LEVEL)
                       .advance(-1)
                       .advance_wrap(1));

        // Check basic properties of advance_wrap().
        // "advancing 7 steps around cube should end up one past start.",
        assert_eq!(CellID::from_face(1),
                   CellID::from_face(0).child_begin_at_level(0).advance_wrap(7));
        // "twice around should end up where we started",
        assert_eq!(CellID::from_face(0).child_begin_at_level(0),
                   CellID::from_face(0)
                       .child_begin_at_level(0)
                       .advance_wrap(12));
        // "backwards once around plus one step should be one before we started",
        assert_eq!(CellID::from_face(4), CellID::from_face(5).advance_wrap(-7));
        // "wrapping even multiple of times around should end where we started",
        assert_eq!(CellID::from_face(0).child_begin_at_level(0),
                   CellID::from_face(0)
                       .child_begin_at_level(0)
                       .advance_wrap(-12000000));
        // "wrapping combination of even times around should end where it started",
        assert_eq!(CellID::from_face(0)
                       .child_begin_at_level(5)
                       .advance_wrap(6644),
                   CellID::from_face(0)
                       .child_begin_at_level(5)
                       .advance_wrap(-11788));
        // "moving 256 should advance us one cell at max level",
        assert_eq!(id.next().child_begin_at_level(MAX_LEVEL),
                   id.child_begin_at_level(MAX_LEVEL).advance_wrap(256));
        // "wrapping by 4 times cells per face should advance 4 faces",
        assert_eq!(CellID::from_face_pos_level(1, 0, MAX_LEVEL),
                   CellID::from_face_pos_level(5, 0, MAX_LEVEL).advance_wrap(2 << (2 * MAX_LEVEL)));
    }

    #[test]
    fn test_cellid_advance() {
        assert_eq!(CellID::from_face(0).child_begin_at_level(0).advance(7),
                   CellID::from_face(5).child_end_at_level(0));
        assert_eq!(CellID::from_face(0).child_begin_at_level(0).advance(12),
                   CellID::from_face(5).child_end_at_level(0));
        assert_eq!(CellID::from_face(5).child_end_at_level(0).advance(-7),
                   CellID::from_face(0).child_begin_at_level(0));
        assert_eq!(CellID::from_face(5)
                       .child_end_at_level(0)
                       .advance(-12000000),
                   CellID::from_face(0).child_begin_at_level(0));
        assert_eq!(CellID::from_face(0).child_begin_at_level(5).advance(500),
                   CellID::from_face(5)
                       .child_end_at_level(5)
                       .advance(500 - (6 << (2 * 5))));
        assert_eq!(CellID::from_face_pos_level(3, 0x12345678, MAX_LEVEL - 4)
                       .child_begin_at_level(MAX_LEVEL)
                       .advance(256),
                   CellID::from_face_pos_level(3, 0x12345678, MAX_LEVEL - 4)
                       .next()
                       .child_begin_at_level(MAX_LEVEL));
        assert_eq!(CellID::from_face_pos_level(1, 0, MAX_LEVEL).advance(4 << (2 * MAX_LEVEL)),
                   CellID::from_face_pos_level(5, 0, MAX_LEVEL));
    }

    #[test]
    fn test_cellid_face_siti() {
        let id = CellID::from_face_pos_level(3, 0x12345678, MAX_LEVEL);

        // Check that the (si, ti) coordinates of the center end in a
        // 1 followed by (30 - level) 0's.
        for level in 0..(MAX_LEVEL + 1) {
            let l = MAX_LEVEL - level;
            let want = 1 << level;
            let mask = (1 << (level + 1)) - 1;

            let (_, si, ti) = id.parent(l).face_siti();
            assert_eq!(want, (si as u64) & mask);
            assert_eq!(want, (ti as u64) & mask);
        }
    }
}


/*
func TestCellIDContinuity(t *testing.T) {
	const maxWalkLevel = 8
	const cellSize = 1.0 / (1 << maxWalkLevel)

	// Make sure that sequentially increasing cell ids form a continuous
	// path over the surface of the sphere, i.e. there are no
	// discontinuous jumps from one region to another.

	maxDist := MaxWidthMetric.Value(maxWalkLevel)
	end := CellID::from_face(5).child_end_at_level(maxWalkLevel)
	id := CellID::from_face(0).child_begin_at_level(maxWalkLevel)

	for ; id != end; id = id.next() {

		if got := id.rawPoint().Angle(id.next_wrap().rawPoint()); float64(got) > maxDist {
			t.Errorf("%v.rawPoint().Angle(%v.next_wrap().rawPoint()) = %v > %v", id, id, got, maxDist)
		}
		if id.next_wrap() != id.advance_wrap(1) {
			t.Errorf("%v.next_wrap() != %v.advance_wrap(1) %v != %v)", id, id, id.next_wrap(), id.advance_wrap(1))
		}
		if id != id.next_wrap().advance_wrap(-1) {
			t.Errorf("%v.next_wrap().advance_wrap(-1) = %v want %v)", id, id.next_wrap().advance_wrap(-1), id)
		}

		// Check that the rawPoint() returns the center of each cell
		// in (s,t) coordinates.
		_, u, v := xyzToFaceUV(id.rawPoint())
		if !float64Eq(math.Remainder(uvToST(u), 0.5*cellSize), 0.0) {
			t.Errorf("uvToST(%v) = %v, want %v", u, uvToST(u), 0.5*cellSize)
		}
		if !float64Eq(math.Remainder(uvToST(v), 0.5*cellSize), 0.0) {
			t.Errorf("uvToST(%v) = %v, want %v", v, uvToST(v), 0.5*cellSize)
		}
	}
}

// sampleBoundary returns a random point on the boundary of the given rectangle.
func sampleBoundary(rect r2.Rect) (u, v float64) {
	if oneIn(2) {
		v = randomUniformFloat64(rect.Y.Lo, rect.Y.Hi)
		if oneIn(2) {
			u = rect.X.Lo
		} else {
			u = rect.X.Hi
		}
	} else {
		u = randomUniformFloat64(rect.X.Lo, rect.X.Hi)
		if oneIn(2) {
			v = rect.Y.Lo
		} else {
			v = rect.Y.Hi
		}
	}
	return u, v
}

// projectToBoundary returns the closest point to uv on the boundary of rect.
func projectToBoundary(u, v float64, rect r2.Rect) r2.Point {
	du0 := math.Abs(u - rect.X.Lo)
	du1 := math.Abs(u - rect.X.Hi)
	dv0 := math.Abs(v - rect.Y.Lo)
	dv1 := math.Abs(v - rect.Y.Hi)

	dmin := math.Min(math.Min(du0, du1), math.Min(dv0, dv1))
	if du0 == dmin {
		return r2.Point{rect.X.Lo, rect.Y.ClampPoint(v)}
	}
	if du1 == dmin {
		return r2.Point{rect.X.Hi, rect.Y.ClampPoint(v)}
	}
	if dv0 == dmin {
		return r2.Point{rect.X.ClampPoint(u), rect.Y.Lo}
	}

	return r2.Point{rect.X.ClampPoint(u), rect.Y.Hi}
}

func TestCellIDExpandedByDistanceUV(t *testing.T) {
	const maxDistDegrees = 10
	for i := 0; i < 1000; i++ {
		id := randomCellID()
		distance := s1.Degree * s1.Angle(randomUniformFloat64(-maxDistDegrees, maxDistDegrees))

		bound := id.boundUV()
		expanded := expandedByDistanceUV(bound, distance)
		for iter := 0; iter < 10; iter++ {
			// Choose a point on the boundary of the rectangle.
			face := randomUniformInt(6)
			centerU, centerV := sampleBoundary(bound)
			center := Point{faceUVToXYZ(face, centerU, centerV).Normalize()}

			// Now sample a point from a disc of radius (2 * distance).
			p := samplePointFromCap(CapFromCenterHeight(center, 2*math.Abs(float64(distance))))

			// Find the closest point on the boundary to the sampled point.
			u, v, ok := faceXYZToUV(face, p)
			if !ok {
				continue
			}

			uv := r2.Point{u, v}
			closestUV := projectToBoundary(u, v, bound)
			closest := faceUVToXYZ(face, closestUV.X, closestUV.Y).Normalize()
			actualDist := p.Distance(Point{closest})

			if distance >= 0 {
				// expanded should contain all points in the original bound,
				// and also all points within distance of the boundary.
				if bound.ContainsPoint(uv) || actualDist < distance {
					if !expanded.ContainsPoint(uv) {
						t.Errorf("expandedByDistanceUV(%v, %v).ContainsPoint(%v) = false, want true", bound, distance, uv)
					}
				}
			} else {
				// expanded should not contain any points within distance
				// of the original boundary.
				if actualDist < -distance {
					if expanded.ContainsPoint(uv) {
						t.Errorf("negatively expandedByDistanceUV(%v, %v).ContainsPoint(%v) = true, want false", bound, distance, uv)
					}
				}
			}
		}
	}
}

func TestCellIDMaxTile(t *testing.T) {
	// This method is also tested more thoroughly in s2cellunion_test.
	for iter := 0; iter < 1000; iter++ {
		id := randomCellIDForLevel(10)

		// Check that limit is returned for tiles at or beyond limit.
		if got, want := id, id.MaxTile(id); got != want {
			t.Errorf("%v.MaxTile(%v) = %v, want %v", id, id, got, want)
		}
		if got, want := id, id.Children()[0].MaxTile(id); got != want {
			t.Errorf("%v.Children()[0].MaxTile(%v) = %v, want %v", id, id, got, want)
		}
		if got, want := id, id.Children()[1].MaxTile(id); got != want {
			t.Errorf("%v.Children()[1].MaxTile(%v) = %v, want %v", id, id, got, want)
		}
		if got, want := id, id.next().MaxTile(id); got != want {
			t.Errorf("%v.next().MaxTile(%v) = %v, want %v", id, id, got, want)
		}
		if got, want := id.Children()[0], id.MaxTile(id.Children()[0]); got != want {
			t.Errorf("%v.MaxTile(%v.Children()[0] = %v, want %v", id, id, got, want)
		}

		// Check that the tile size is increased when possible.
		if got, want := id, id.Children()[0].MaxTile(id.next()); got != want {
			t.Errorf("%v.Children()[0].MaxTile(%v.next()) = %v, want %v", id, id, got, want)
		}

		if got, want := id, id.Children()[0].MaxTile(id.next().Children()[0]); got != want {
			t.Errorf("%v.Children()[0].MaxTile(%v.next()) = %v, want %v", id, id, got, want)
		}

		if got, want := id, id.Children()[0].MaxTile(id.next().Children()[1].Children()[0]); got != want {
			t.Errorf("%v.Children()[0].MaxTile(%v.next().Children()[1].Children()[0] = %v, want %v", id, id, got, want)
		}

		if got, want := id, id.Children()[0].Children()[0].MaxTile(id.next()); got != want {
			t.Errorf("%v.Children()[0].Children()[0].MaxTile(%v.next()) = %v, want %v", id, id, got, want)
		}

		if got, want := id, id.Children()[0].Children()[0].Children()[0].MaxTile(id.next()); got != want {
			t.Errorf("%v.Children()[0].Children()[0].Children()[0].MaxTile(%v.next()) = %v, want %v", id, id, got, want)
		}

		// Check that the tile size is decreased when necessary.
		if got, want := id.Children()[0], id.MaxTile(id.Children()[0].next()); got != want {
			t.Errorf("%v.Children()[0], id.MaxTile(%v.Children()[0].next()) = %v, want %v", id, id, got, want)
		}

		if got, want := id.Children()[0], id.MaxTile(id.Children()[0].next().Children()[0]); got != want {
			t.Errorf("%v.Children()[0], id.MaxTile(%v.Children()[0].next().Children()[0]) = %v, want %v", id, id, got, want)
		}

		if got, want := id.Children()[0], id.MaxTile(id.Children()[0].next().Children()[1]); got != want {
			t.Errorf("%v.Children()[0], id.MaxTile(%v.Children()[0].next().Children()[1]) = %v, want %v", id, id, got, want)
		}

		if got, want := id.Children()[0].Children()[0], id.MaxTile(id.Children()[0].Children()[0].next()); got != want {
            t.Errorf("%v.Children()[0].Children()[0],
                     id.MaxTile(%v.Children()[0].Children()[0].next()) = %v, want %v", id, id, got,
                     want)
		}

		if got, want := id.Children()[0].Children()[0].Children()[0],
			id.MaxTile(id.Children()[0].Children()[0].Children()[0].next()); got != want {
			t.Errorf("%v.MaxTile(%v.Children()[0].Children()[0].Children()[0].next()) = %v, want %v", id, id, got, want)
		}

		// Check that the tile size is otherwise unchanged.
		if got, want := id, id.MaxTile(id.next()); got != want {
			t.Errorf("%v.MaxTile(%v.next()) = %v, want %v", id, id, got, want)
		}

		if got, want := id, id.MaxTile(id.next().Children()[0]); got != want {
			t.Errorf("%v.MaxTile(%v.next().Children()[0]) = %v, want %v", id, id, got, want)
		}

		if got, want := id, id.MaxTile(id.next().Children()[1].Children()[0]); got != want {
			t.Errorf("%v.MaxTile(%v.next().Children()[1].Children()[0]) = %v, want %v", id, id, got, want)
		}
	}
}

func TestCellIDCenterFaceSiTi(t *testing.T) {
	// Check that the (si, ti) coordinates of the center end in a
	// 1 followed by (30 - level) 0s.

	id := CellID::from_face_pos_level(3, 0x12345678, MAX_LEVEL)

	tests := []struct {
		id          CellID
		levelOffset uint
	}{
		// Leaf level, 30.
		{id, 0},
		// Level 29.
		{id.parent(MAX_LEVEL - 1), 1},
		// Level 28.
		{id.parent(MAX_LEVEL - 2), 2},
		// Level 20.
		{id.parent(MAX_LEVEL - 10), 10},
		// Level 10.
		{id.parent(MAX_LEVEL - 20), 20},
		// Level 0.
		{id.parent(0), MAX_LEVEL},
	}

	for _, test := range tests {
		_, si, ti := test.id.centerFaceSiTi()
		want := 1 << test.levelOffset
		mask := (1 << (test.levelOffset + 1)) - 1
		if want != si&mask {
			t.Errorf("Level Offset %d. %b != %b", test.levelOffset, want, si&mask)
		}
		if want != ti&mask {
			t.Errorf("Level Offset: %d. %b != %b", test.levelOffset, want, ti&mask)
		}
	}
}

// TODO(roberts): Remaining tests to convert.
// Coverage
// TraversalOrder
*/
