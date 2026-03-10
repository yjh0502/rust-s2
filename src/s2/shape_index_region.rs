// Copyright 2023 Google Inc. All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use crate::s2::cap::Cap;
use crate::s2::cellid::CellID;
use crate::s2::cellunion::CellUnion;
use crate::s2::point::Point;
use crate::s2::rect::Rect;
use crate::s2::region::Region;
use crate::s2::shape_index::{ShapeIndex, ShapeIndexIterator};

/// A wrapper used for ContainsPointQuery. This will need to be properly implemented.
/// This is a temporary implementation to make the code compile.
#[derive(Clone)]
pub struct ContainsPointQuery<'a> {
    index: &'a ShapeIndex,
    vertex_model: VertexModel,
}

impl<'a> ContainsPointQuery<'a> {
    pub fn new(index: &'a ShapeIndex, vertex_model: VertexModel) -> Self {
        ContainsPointQuery {
            index,
            vertex_model,
        }
    }

    pub fn contains(&self, _point: Point) -> bool {
        // This is a placeholder implementation
        // TODO: Implement properly
        false
    }
}

/// Defines various ways of handling vertices for polygon containment tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexModel {
    /// All vertices are contained within the shape.
    Open,
    /// No vertices are contained within the shape.
    Closed,
    /// Vertices that are contained within the shape are defined using symbolic perturbations.
    SemiOpen,
}

/// ShapeIndexRegion wraps a ShapeIndex and implements the Region interface.
/// This allows RegionCoverer to work with ShapeIndexes as well as being
/// able to be used by some of the Query types.
#[derive(Clone)]
pub struct ShapeIndexRegion<'a> {
    pub index: &'a ShapeIndex,
    pub contains_query: ContainsPointQuery<'a>,
    pub iter: ShapeIndexIterator<'a>,
}

// TODO: Uncomment once implementation is complete.
// Enforce Region interface satisfaction similar to other types that implement Region.
// impl<'a> Region for ShapeIndexRegion<'a> {}

impl<'a> Region for ShapeIndexRegion<'a> {
    /// cap_bound returns a bounding spherical cap for this collection of geometry.
    /// This is not guaranteed to be exact.
    fn cap_bound(&self) -> Cap {
        let cu = CellUnion(self.cell_union_bound());
        cu.cap_bound()
    }

    /// rect_bound returns a bounding rectangle for this collection of geometry.
    /// The bounds are not guaranteed to be tight.
    fn rect_bound(&self) -> Rect {
        let cu = CellUnion(self.cell_union_bound());
        cu.rect_bound()
    }

    // TODO: Implement additional Region interface methods
    // fn contains_cell(&self, target: &Cell) -> bool { ... }
    // fn intersects_cell(&self, target: &Cell) -> bool { ... }
}

impl<'a> ShapeIndexRegion<'a> {
    /// cell_union_bound returns the bounding CellUnion for this collection of geometry.
    /// This method currently returns at most 4 cells, unless the index spans
    /// multiple faces in which case it may return up to 6 cells.
    pub fn cell_union_bound(&self) -> Vec<CellID> {
        // We find the range of Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  There are
        // two cases:
        //
        //  - If the index intersects two or more faces, then for each intersected
        //    face we add one cell to the covering.  Rather than adding the entire
        //    face, instead we add the smallest Cell that covers the ShapeIndex
        //    cells within that face.
        //
        //  - If the index intersects only one face, then we first find the smallest
        //    cell S that contains the index cells (just like the case above).
        //    However rather than using the cell S itself, instead we repeat this
        //    process for each of its child cells.  In other words, for each
        //    child cell C we add the smallest Cell C' that covers the index cells
        //    within C.  This extra step is relatively cheap and produces much
        //    tighter coverings when the ShapeIndex consists of a small region
        //    near the center of a large Cell.
        let mut cell_ids = Vec::new();

        // Create a copy of the iterator to avoid modifying the original
        let mut iterator = self.iter.clone();

        // Find the last CellID in the index.
        iterator.end();
        if !iterator.prev() {
            return cell_ids; // Empty index.
        }

        let last_index_id = iterator.cell_id();
        iterator.begin();

        if iterator.cell_id() != last_index_id {
            // The index has at least two cells. Choose a CellID level such that
            // the entire index can be spanned with at most 6 cells (if the index
            // spans multiple faces) or 4 cells (it the index spans a single face).
            let level_result = iterator.cell_id().common_ancestor_level(&last_index_id);

            let mut level = match level_result {
                Some(l) => l,
                None => 0, // No common ancestor, equivalent to level -1 in Go
            };

            // In Go, level is incremented after setting to -1 when no common ancestor
            // In Rust, we've already set it to 0, so we increment to 1
            level += 1;

            // For each cell C at the chosen level, we compute the smallest Cell
            // that covers the ShapeIndex cells within C.
            let last_id = last_index_id.parent(level);

            let mut id = iterator.cell_id().parent(level);
            while id != last_id {
                // If the cell C does not contain any index cells, then skip it.
                if id.range_max() < iterator.cell_id() {
                    id = id.next();
                    continue;
                }

                // Find the range of index cells contained by C and then shrink C so
                // that it just covers those cells.
                let first = iterator.cell_id();
                iterator.seek(id.range_max().next());
                iterator.prev();

                cell_ids = self.cover_range(first, iterator.cell_id(), cell_ids);
                iterator.next();

                id = id.next();
            }
        }

        self.cover_range(iterator.cell_id(), last_index_id, cell_ids)
    }

    /// cover_range computes the smallest CellID that covers the Cell range (first, last)
    /// and returns the updated slice.
    ///
    /// This requires first and last have a common ancestor.
    fn cover_range(&self, first: CellID, last: CellID, mut cell_ids: Vec<CellID>) -> Vec<CellID> {
        // The range consists of a single index cell.
        if first == last {
            cell_ids.push(first);
            return cell_ids;
        }

        // Add the lowest common ancestor of the given range.
        match first.common_ancestor_level(&last) {
            Some(level) => {
                cell_ids.push(first.parent(level));
                cell_ids
            }
            None => {
                // No common ancestor
                cell_ids.push(CellID(0));
                cell_ids
            }
        }
    }

    // TODO: Implement remaining methods
    /*
    pub fn contains_cell(&self, target: &Cell) -> bool {
        // ...
    }

    pub fn intersects_cell(&self, target: &Cell) -> bool {
        // ...
    }

    pub fn contains_point(&self, p: Point) -> bool {
        // ...
    }

    fn contains(&self, id: CellID, clipped: &ClippedShape, p: Point) -> bool {
        // ...
    }

    fn any_edge_intersects(&self, clipped: &ClippedShape, target: &Cell) -> bool {
        // ...
    }
    */
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::consts::DBL_EPSILON;
    use crate::r2;
    use crate::r2::rect::Rect as R2Rect;
    use std::convert::TryInto;
    // use crate::r2::point::Point as R2Point;
    use crate::s2::cell::Cell;
    use crate::s2::cellid::CellID;
    use crate::s2::lax_loop::LaxLoop;

    use crate::s2::stuv::face_uv_to_xyz;
    use crate::shape::ShapeType;
    use crate::shape_index_region::ShapeIndex;

    // Set padding to at least twice the maximum error for reliable results.
    // TODO: Use the constants from edgeutil.rs once fully ported
    const FACE_CLIP_ERROR_UV_COORD: f64 = 9.0 * (1.0 / std::f64::consts::SQRT_2) * DBL_EPSILON;
    const INTERSECTS_RECT_ERROR_UV_DIST: f64 = 3.0 * DBL_EPSILON;
    const SHAPE_INDEX_CELL_PADDING: f64 =
        2.0 * (FACE_CLIP_ERROR_UV_COORD + INTERSECTS_RECT_ERROR_UV_DIST);

    // Helper function to convert from face, i, j coordinates to uv-coordinates
    fn ij_level_to_bound_uv(i: u64, j: u64, level: u64) -> R2Rect {
        // This is a simplified implementation
        // In a complete implementation, this would convert i,j coordinates to UV bounds
        let size = 1.0 / (1u64 << level) as f64;
        let u = (i as f64) * size;
        let v = (j as f64) * size;

        R2Rect::from_points(&[
            r2::point::Point::new(u, v),
            r2::point::Point::new(u + size, v),
            r2::point::Point::new(u + size, v + size),
            r2::point::Point::new(u, v + size),
        ])
    }

    // Pad a cell with the given UV padding
    fn pad_cell(id: CellID, padding_uv: f64) -> ShapeType {
        let (face, i, j, _) = id.face_ij_orientation();

        // In the Go code, this is uv := ijLevelToBoundUV(i, j, id.Level()).ExpandedByMargin(paddingUV)
        let uv = ij_level_to_bound_uv(i.try_into().unwrap(), j as u64, id.level())
            .expanded_by_margin(padding_uv);

        let mut vertices = Vec::with_capacity(4);
        for vertex in uv.vertices() {
            let xyz = face_uv_to_xyz(face, vertex.x, vertex.y).normalize();
            vertices.push(xyz.into());
        }

        LaxLoop::from_points(vertices).into()
    }

    // Helper function to sort CellIDs for comparison
    fn sort_cell_ids(ids: &mut Vec<CellID>) {
        ids.sort();
    }

    #[test]
    fn test_shape_index_region_cap_bound() {
        // In Go: id := CellIDFromString("3/0123012301230123012301230123")
        let id = CellID::from_string("3/0123012301230123012301230123");

        // Add a polygon that is slightly smaller than the cell being tested
        let mut index = ShapeIndex::new();
        index.add(&pad_cell(id, -SHAPE_INDEX_CELL_PADDING));

        let cell_bound = Cell::from(&id).cap_bound();
        let index_bound = index.region().cap_bound();

        assert!(
            index_bound.contains(&cell_bound),
            "Expected index_bound {:?} to contain cell_bound {:?}",
            index_bound,
            cell_bound
        );

        // Note that CellUnion.CapBound returns a slightly larger bound than
        // Cell.CapBound even when the cell union consists of a single CellID.
        let got = index_bound.radius();
        let want = 1.00001 * cell_bound.radius();
        assert!(
            got <= want,
            "index.CapBound.Radius() = {:?}, want <= {:?}",
            got,
            want
        );
    }

    #[test]
    fn test_shape_index_region_rect_bound() {
        // In Go: id := CellIDFromString("3/0123012301230123012301230123")
        let id = CellID::from_string("3/0123012301230123012301230123");

        // Add a polygon that is slightly smaller than the cell being tested
        let mut index = ShapeIndex::new();
        index.add(&pad_cell(id, -SHAPE_INDEX_CELL_PADDING));

        let cell_bound = Cell::from(&id).rect_bound();
        let index_bound = index.region().rect_bound();

        assert_eq!(
            index_bound, cell_bound,
            "index.RectBound() = {:?}, want {:?}",
            index_bound, cell_bound
        );
    }

    #[test]
    fn test_shape_index_region_cell_union_bound_multiple_faces() {
        let have = vec![
            CellID::from_string("3/00123"),
            CellID::from_string("2/11200013"),
        ];

        let mut index = ShapeIndex::new();
        for id in &have {
            index.add(&pad_cell(*id, -SHAPE_INDEX_CELL_PADDING));
        }

        let mut got = index.region().cell_union_bound();

        // Sort both for comparison
        let mut have_copy = have.clone();
        sort_cell_ids(&mut have_copy);
        sort_cell_ids(&mut got);

        assert_eq!(
            CellUnion(got.clone()),
            CellUnion(have_copy.clone()),
            "index.CellUnionBound() = {:?}, want {:?}",
            got,
            have_copy
        );
    }

    #[test]
    fn test_shape_index_region_cell_union_bound_one_face() {
        // This tests consists of 3 pairs of CellIDs.  Each pair is located within
        // one of the children of face 5, namely the cells 5/0, 5/1, and 5/3.
        // We expect CellUnionBound to compute the smallest cell that bounds the
        // pair on each face.
        let have = vec![
            CellID::from_string("5/010"),
            CellID::from_string("5/0211030"),
            CellID::from_string("5/110230123"),
            CellID::from_string("5/11023021133"),
            CellID::from_string("5/311020003003030303"),
            CellID::from_string("5/311020023"),
        ];

        let want = vec![
            CellID::from_string("5/0"),
            CellID::from_string("5/110230"),
            CellID::from_string("5/3110200"),
        ];

        let mut index = ShapeIndex::new();
        for id in &have {
            // Add each shape 3 times to ensure that the ShapeIndex subdivides.
            index.add(&pad_cell(*id, -SHAPE_INDEX_CELL_PADDING));
            index.add(&pad_cell(*id, -SHAPE_INDEX_CELL_PADDING));
            index.add(&pad_cell(*id, -SHAPE_INDEX_CELL_PADDING));
        }

        let mut have_copy = have.clone();
        sort_cell_ids(&mut have_copy);

        let got = index.region().cell_union_bound();

        assert_eq!(
            CellUnion(got.clone()),
            CellUnion(want.clone()),
            "index.CellUnionBound() = {:?}, want {:?}",
            got,
            want
        );
    }

    // TODO: Implement remaining tests
    // fn test_shape_index_region_contains_cell_multiple_shapes() { ... }
    // fn test_shape_index_region_intersects_shrunken_cell() { ... }
    // fn test_shape_index_region_intersects_exact_cell() { ... }
}
