/*
Copyright 2015 Google Inc. All rights reserved.
Copyright 2023-2024 Jhyun Yu and contributors. All rights reserved.

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
use super::region::Region;
use crate::r#loop::{BoundaryCondition, Loop};
use crate::r3::vector::Vector as R3Vector;
use crate::s2::cap::Cap;
use crate::s2::cell::Cell;
use crate::s2::edge_clipping::{
    clip_to_padded_face, edge_intersects_rect, INTERSECT_RECT_ERROR_UV_DIST,
};
use crate::s2::edge_crosser::EdgeCrosser;
use crate::s2::point::Point;
use crate::s2::rect::Rect;
use crate::s2::rect_bounder::expand_for_subregions;
use crate::s2::shape::{Chain, ChainPosition, Edge, ReferencePoint, Shape};
use crate::s2::shape_index::{ShapeIndex, ShapeIndexIterator};
use crate::shape::ShapeType;
use crate::shape_index::CellRelation::{Disjoint, Indexed, Subdivided};
use crate::shape_index::FACE_CLIP_ERROR_UV_COORD;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::TryInto;
use std::hash::{Hash, Hasher};

/// Polygon represents a sequence of zero or more loops; recall that the
/// interior of a loop is defined to be its left-hand side (see Loop).
///
/// When the polygon is initialized, the given loops are automatically converted
/// into a canonical form consisting of "shells" and "holes". Shells and holes
/// are both oriented CCW, and are nested hierarchically. The loops are
/// reordered to correspond to a pre-order traversal of the nesting hierarchy.
///
/// Polygons may represent any region of the sphere with a polygonal boundary,
/// including the entire sphere (known as the "full" polygon). The full polygon
/// consists of a single full loop (see Loop), whereas the empty polygon has no
/// loops at all.
///
/// Use full_polygon() to construct a full polygon. The default constructor creates an empty polygon.
///
/// Polygons have the following restrictions:
///
///   - Loops may not cross, i.e. the boundary of a loop may not intersect
///     both the interior and exterior of any other loop.
///
///   - Loops may not share edges, i.e. if a loop contains an edge AB, then
///     no other loop may contain AB or BA.
///
///   - Loops may share vertices, however no vertex may appear twice in a
///     single loop (see Loop).
///
///   - No loop may be empty. The full loop may appear only in the full polygon.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Polygon {
    loops: Vec<Loop>,

    // index is a spatial index of all the polygon loops.
    #[cfg_attr(feature = "serde", serde(skip))]
    index: ShapeIndex,

    // has_holes tracks if this polygon has at least one hole.
    has_holes: bool,

    // num_vertices keeps the running total of all of the vertices of the contained loops.
    num_vertices: usize,

    // num_edges tracks the total number of edges in all the loops in this polygon.
    num_edges: usize,

    // bound is a conservative bound on all points contained by this loop.
    // If l.contains_point(P), then l.bound.contains_point(P).
    bound: Rect,

    // Since bound is not exact, it is possible that a loop A contains
    // another loop B whose bounds are slightly larger. subregion_bound
    // has been expanded sufficiently to account for this error, i.e.
    // if A.contains(B), then A.subregion_bound.contains(B.bound).
    subregion_bound: Rect,

    // A vector where element i is the cumulative number of edges in the
    // preceding loops in the polygon. This field is used for polygons that
    // have a large number of loops, and may be empty for polygons with few loops.
    cumulative_edges: Vec<usize>,
}

// Loop relation is a map of a loop to its immediate children with respect to nesting.
// It is used to determine which loops are shells and which are holes.
type LoopMap<'a> = HashMap<&'a Loop, Vec<&'a Loop>>;

impl Default for Polygon {
    fn default() -> Self {
        Self {
            loops: Vec::new(),
            index: ShapeIndex::new(),
            has_holes: false,
            num_vertices: 0,
            num_edges: 0,
            bound: Rect::empty(),
            subregion_bound: Rect::empty(),
            cumulative_edges: Vec::new(),
        }
    }
}

impl Hash for Polygon {
    fn hash<H: Hasher>(&self, _state: &mut H) {
        todo!()
    }
}

impl PartialEq for Polygon {
    fn eq(&self, _other: &Self) -> bool {
        todo!()
    }
}

impl Eq for Polygon {}

// Implementation of basic Polygon methods
impl Polygon {
    /// Creates a new empty polygon
    pub fn new() -> Self {
        Self::default()
    }

    /// from_loops constructs a polygon from the given set of loops. The polygon
    /// interior consists of the points contained by an odd number of loops. (Recall
    /// that a loop contains the set of points on its left-hand side.)
    ///
    /// This method determines the loop nesting hierarchy and assigns every loop a
    /// depth. Shells have even depths, and holes have odd depths.
    ///
    /// Note: The given set of loops are reordered by this method so that the hierarchy
    /// can be traversed using parent, last_descendant and the loops depths.
    pub fn from_loops(loops: Vec<Loop>) -> Self {
        let mut p = Self::default();

        // Empty polygons do not contain any loops, even the Empty loop.
        if loops.len() == 1 && loops[0].is_empty() {
            p.init_loop_properties();
            return p;
        }

        p.loops = loops;
        p.init_nested();
        p
    }

    /// from_oriented_loops returns a Polygon from the given set of loops,
    /// like from_loops. It expects loops to be oriented such that the polygon
    /// interior is on the left-hand side of all loops. This implies that shells
    /// and holes should have opposite orientations in the input to this method.
    /// (During initialization, loops representing holes will automatically be
    /// inverted.)
    pub fn from_oriented_loops(mut loops: Vec<Loop>) -> Self {
        // Remember which of the given loops contain the origin
        let mut contained_origin = HashMap::new();
        for loop_ref in &loops {
            contained_origin.insert(loop_ref as *const Loop, loop_ref.contains_origin());
        }

        for loop_ref in &mut loops {
            let angle = loop_ref.turning_angle();
            if angle.abs() > loop_ref.turning_angle_max_error() {
                // Normalize the loop.
                if angle < 0.0 {
                    loop_ref.invert();
                }
            } else {
                // Ensure that the loop does not contain the origin.
                if loop_ref.contains_origin() {
                    loop_ref.invert();
                }
            }
        }

        let mut p = Self::from_loops(loops);

        if p.num_loops() > 0 {
            let mut origin_loop = &p.loops[0];
            let mut polygon_contains_origin = false;

            for l in &p.loops {
                if l.contains_origin() {
                    polygon_contains_origin = !polygon_contains_origin;
                    origin_loop = l;
                }
            }

            if contained_origin.get(&(origin_loop as *const Loop)) != Some(&polygon_contains_origin)
            {
                p.invert();
            }
        }

        p
    }

    /// from_cell returns a Polygon from a single loop created from the given Cell.
    pub fn from_cell(cell: &Cell) -> Self {
        Self::from_loops(vec![Loop::from_cell(cell)])
    }

    /// init_nested takes the set of loops in this polygon and performs the nesting
    /// computations to set the proper nesting and parent/child relationships.
    fn init_nested(&mut self) {
        if self.loops.len() == 1 {
            self.init_one_loop();
            return;
        }

        // Take ownership of loops to avoid borrow issues
        let loops = std::mem::take(&mut self.loops);

        // Setup loop map for ordering
        let mut lm = LoopMap::new();

        // Create a root entry in the map to represent "no parent"
        // Insert loops into loop map for ordering
        Self::build_loop_hierarchy(&mut lm, &loops);

        // Reorder the loops in depth-first traversal order
        self.init_loops(&lm);
        self.init_loop_properties();
    }

    /// build_loop_hierarchy constructs the parent-child relationship hierarchy
    /// for all loops in the polygon.
    fn build_loop_hierarchy<'a>(lm: &mut LoopMap<'a>, loops: &'a [Loop]) {
        // First pass: add all loops to the map with empty children lists
        for l in loops {
            lm.insert(l, Vec::new());
        }

        // Second pass: for each loop, find its parent and update the hierarchy
        for new_loop in loops {
            Self::find_parent_and_update(lm, new_loop);
        }
    }

    /// find_parent_and_update finds the appropriate parent for a loop and updates
    /// the hierarchy accordingly.
    fn find_parent_and_update<'a>(lm: &mut LoopMap<'a>, new_loop: &'a Loop) {
        // Find the innermost loop that contains new_loop
        let mut parent = None;

        // First, check which loops could contain this one
        for &potential_parent in lm.keys() {
            // Skip the loop itself
            if std::ptr::eq(potential_parent, new_loop) {
                continue;
            }

            // If this loop contains the new loop and is nested deeper than our current parent
            if potential_parent.contains_nested(new_loop) {
                match parent {
                    None => parent = Some(potential_parent),
                    Some(current_parent) => {
                        // If this potential parent is contained by our current parent,
                        // it's a better (more specific) parent
                        if current_parent.contains_nested(potential_parent) {
                            parent = Some(potential_parent);
                        }
                    }
                }
            }
        }

        // Get the parent loop or use None if there's no parent
        if let Some(parent_loop) = parent {
            // Check if any existing children of the parent should now be children of new_loop
            let mut current_children = lm.get(parent_loop).cloned().unwrap_or_default();
            let mut new_children = lm.get(new_loop).cloned().unwrap_or_default();

            // Move appropriate children
            let mut i = 0;
            while i < current_children.len() {
                let child = current_children[i];
                // Skip the new loop itself
                if std::ptr::eq(child, new_loop) {
                    i += 1;
                    continue;
                }

                if new_loop.contains_nested(child) {
                    new_children.push(child);
                    current_children.remove(i);
                } else {
                    i += 1;
                }
            }

            // Add the new loop as a child of the parent
            current_children.push(new_loop);

            // Update the loop map
            lm.insert(parent_loop, current_children);
            lm.insert(new_loop, new_children);
        } else {
            // This is a top-level loop with no parent
            // We don't need to modify its children, they were already set in the initial pass
        }
    }

    /// init_loops walks the mapping of loops to all of their children, and adds them in
    /// order into to the polygons set of loops.
    fn init_loops(&mut self, lm: &LoopMap) {
        let mut stack: Vec<&Loop> = vec![];
        let mut depth = -1;

        while !stack.is_empty() {
            let loop_target = stack.pop();

            if let Some(loop_target) = loop_target {
                depth = (*loop_target).depth;
                self.loops.push((*loop_target).clone());
            }

            // TODO: This is a odd enough pattern in go already (getting a potentially null pointer)
            //       perhaps it is better to refactor into better rust idiom once the tests are good.
            let children = lm.get(&loop_target.unwrap()).cloned().unwrap_or_default();

            // Add children in reverse order so they're processed in forward order
            for child in children.iter().rev() {
                let mut loop_clone = unsafe { (**child).clone() };
                loop_clone.depth = depth + 1;
                let clone_ptr = child; // We're using the same pointer
                stack.push(*clone_ptr);
            }
        }
    }

    /// init_one_loop set the properties for a polygon made of a single loop.
    fn init_one_loop(&mut self) {
        self.has_holes = false;
        self.num_vertices = self.loops[0].num_vertices();
        self.bound = self.loops[0].rect_bound();
        self.subregion_bound = expand_for_subregions(&self.bound);
        // Ensure the loops depth is set correctly.
        self.loops[0].depth = 0;

        self.init_edges_and_index();
    }

    /// init_loop_properties sets the properties for polygons with multiple loops.
    fn init_loop_properties(&mut self) {
        self.num_vertices = 0;
        // the loops depths are set by init_nested/initOriented prior to this.
        self.bound = Rect::empty();
        self.has_holes = false;

        for l in &self.loops {
            if l.is_hole() {
                self.has_holes = true;
            } else {
                self.bound = self.bound.union(&l.rect_bound());
            }
            self.num_vertices += l.num_vertices();
        }

        self.subregion_bound = expand_for_subregions(&self.bound);
        self.init_edges_and_index();
    }

    /// init_edges_and_index performs the shape related initializations and adds the final
    /// polygon to the index.
    fn init_edges_and_index(&mut self) {
        self.num_edges = 0;
        self.cumulative_edges = Vec::new();

        if self.is_full() {
            return;
        }

        const MAX_LINEAR_SEARCH_LOOPS: usize = 12; // Based on benchmarks.
        if self.loops.len() > MAX_LINEAR_SEARCH_LOOPS {
            self.cumulative_edges = Vec::with_capacity(self.loops.len());
        }

        for l in &self.loops {
            if !self.cumulative_edges.is_empty() {
                self.cumulative_edges.push(self.num_edges);
            }
            self.num_edges += l.num_vertices();
        }

        self.index = ShapeIndex::new();
        self.index.add(&ShapeType::Polygon(self.clone())); // Note: Need to implement Shape trait for Polygon
    }

    /// full_polygon returns a special "full" polygon.
    pub fn full_polygon() -> Self {
        let mut ret = Self {
            loops: vec![Loop::full()],
            index: ShapeIndex::new(),
            has_holes: false,
            num_vertices: Loop::full().num_vertices(),
            num_edges: 0,
            bound: Rect::full(),
            subregion_bound: Rect::full(),
            cumulative_edges: Vec::new(),
        };
        ret.init_edges_and_index();
        ret
    }

    // Basic accessors

    /// num_loops returns the number of loops in this polygon.
    pub fn num_loops(&self) -> usize {
        self.loops.len()
    }

    /// loops returns a reference to the loops in this polygon.
    pub fn loops(&self) -> &[Loop] {
        &self.loops
    }

    /// loop_at returns the loop at the given index. Note that during initialization,
    /// the given loops are reordered according to a pre-order traversal of the loop
    /// nesting hierarchy. This implies that every loop is immediately followed by
    /// its descendants. This hierarchy can be traversed using the methods parent,
    /// last_descendant, and Loop.depth.
    pub fn loop_at(&self, k: usize) -> &Loop {
        &self.loops[k]
    }

    /// is_empty reports whether this is the special "empty" polygon (consisting of no loops).
    pub fn is_empty(&self) -> bool {
        self.loops.is_empty()
    }

    /// is_full reports whether this is the special "full" polygon (consisting of a
    /// single loop that encompasses the entire sphere).
    pub fn is_full(&self) -> bool {
        self.loops.len() == 1 && self.loops[0].is_full()
    }

    /// invert inverts the polygon (replaces it by its complement).
    pub fn invert(&mut self) {
        // Inverting any one loop will invert the polygon. The best loop to invert
        // is the one whose area is largest, since this yields the smallest area
        // after inversion. The loop with the largest area is always at depth 0.
        // The descendants of this loop all have their depth reduced by 1, while the
        // former siblings of this loop all have their depth increased by 1.

        // The empty and full polygons are handled specially.
        if self.is_empty() {
            *self = Self::full_polygon();
            self.init_loop_properties();
            return;
        }

        if self.is_full() {
            *self = Self::default();
            self.init_loop_properties();
            return;
        }

        // Find the loop whose area is largest (i.e., whose turning angle is
        // smallest), minimizing calls to turning_angle(). In particular, for
        // polygons with a single shell at level 0 there is no need to call
        // turning_angle() at all. (This method is relatively expensive.)
        let mut best = 0;
        const NONE: f64 = 10.0; // Flag that means "not computed yet"
        let mut best_angle = NONE;

        for i in 1..self.num_loops() {
            if self.loop_at(i).depth != 0 {
                continue;
            }

            // We defer computing the turning angle of loop 0 until we discover
            // that the polygon has another top-level shell.
            if best_angle == NONE {
                best_angle = self.loop_at(best).turning_angle();
            }

            let angle = self.loop_at(i).turning_angle();

            // We break ties deterministically in order to avoid having the output
            // depend on the input order of the loops.
            if angle < best_angle
                || (angle == best_angle && compare_loops(self.loop_at(i), self.loop_at(best)) < 0)
            {
                best = i;
                best_angle = angle;
            }
        }

        // Build the new loops vector, starting with the inverted loop.
        let mut best_loop = self.loops[best].clone();
        best_loop.invert();

        let mut new_loops = vec![best_loop];

        // Add the former siblings of this loop as descendants.
        let last_best = self.last_descendant(best);

        for i in 0..self.loops.len() {
            if i < best || i > last_best {
                let mut loop_clone = self.loops[i].clone();
                loop_clone.depth = loop_clone.depth + 1;
                new_loops.push(loop_clone);
            }
        }

        // Add the former children of this loop as siblings.
        for i in 0..self.loops.len() {
            if i > best && i <= last_best {
                let mut loop_clone = self.loops[i].clone();
                loop_clone.depth = loop_clone.depth - 1;
                new_loops.push(loop_clone);
            }
        }

        self.loops = new_loops;
        self.init_loop_properties();
    }

    /// parent returns the index of the parent of loop k.
    /// If the loop does not have a parent, (None) is returned.
    pub fn parent(&self, k: usize) -> Option<usize> {
        // See where we are on the depth hierarchy.
        let depth = self.loops[k].depth;
        if depth == 0 {
            return None;
        }

        // There may be several loops at the same nesting level as us that share a
        // parent loop with us. (Imagine a slice of swiss cheese, of which we are one loop.
        // we don't know how many may be next to us before we get back to our parent loop.)
        // Move up one position from us, and then begin traversing back through the set of loops
        // until we find the one that is our parent or we get to the top of the polygon.
        let mut k = k as isize - 1;
        while k >= 0 && self.loops[k as usize].depth <= depth {
            k -= 1;
        }

        if k < 0 {
            None
        } else {
            Some(k as usize)
        }
    }

    /// last_descendant returns the index of the last loop that is contained within loop k.
    /// If k is negative, it returns the last loop in the polygon.
    /// Note that loops are indexed according to a pre-order traversal of the nesting
    /// hierarchy, so the immediate children of loop k can be found by iterating over
    /// the loops (k+1)..last_descendant(k) and selecting those whose depth is equal
    /// to loop_at(k).depth+1.
    pub fn last_descendant(&self, k: usize) -> usize {
        if k >= self.loops.len() {
            return self.loops.len() - 1;
        }

        let depth = self.loops[k].depth;

        // Find the next loop immediately past us in the set of loops, and then start
        // moving down the list until we either get to the end or find the next loop
        // that is higher up the hierarchy than we are.
        let mut k = k + 1;
        while k < self.loops.len() && self.loops[k].depth > depth {
            k += 1;
        }

        k - 1
    }

    /// any_loop_contains reports whether any loop in this polygon contains the given loop.
    fn any_loop_contains(&self, o: &Loop) -> bool {
        for l in &self.loops {
            if l.contains(o) {
                return true;
            }
        }
        false
    }

    /// any_loop_intersects reports whether any loop in this polygon intersects the given loop.
    fn any_loop_intersects(&self, o: &Loop) -> bool {
        for l in &self.loops {
            if l.intersects(o) {
                return true;
            }
        }
        false
    }

    /// validate checks whether this is a valid polygon,
    /// including checking whether all the loops are themselves valid.
    pub fn validate(&self) -> Result<(), String> {
        for (i, l) in self.loops.iter().enumerate() {
            // Check for loop errors that don't require building a ShapeIndex.
            if let Err(err) = l.find_validation_error_no_index() {
                return Err(format!("loop {}: {}", i, err));
            }

            // Check that no loop is empty, and that the full loop only appears in the
            // full polygon.
            if l.is_empty() {
                return Err(format!("loop {}: empty loops are not allowed", i));
            }

            if l.is_full() && self.loops.len() > 1 {
                return Err(format!("loop {}: full loop appears in non-full polygon", i));
            }
        }

        // TODO: Uncomment the remaining checks when they are completed.

        // Check for loop self-intersections and loop pairs that cross
        // (including duplicate edges and vertices).
        // if findSelfIntersection(p.index) {
        //     return Err("polygon has loop pairs that cross".to_string());
        // }

        // Check whether initOriented detected inconsistent loop orientations.
        // if p.hasInconsistentLoopOrientations {
        //     return Err("inconsistent loop orientations detected".to_string());
        // }

        // Finally, verify the loop nesting hierarchy.
        self.find_loop_nesting_error()
    }

    /// find_loop_nesting_error reports if there is an error in the loop nesting hierarchy.
    fn find_loop_nesting_error(&self) -> Result<(), String> {
        // First check that the loop depths make sense.
        let mut last_depth = -1;
        for (i, l) in self.loops.iter().enumerate() {
            let depth = l.depth;
            if depth < 0 || depth > last_depth + 1 {
                return Err(format!("loop {}: invalid loop depth ({})", i, depth));
            }
            last_depth = depth;
        }

        // Then check that they correspond to the actual loop nesting. This test
        // is quadratic in the number of loops but the cost per iteration is small.
        for i in 0..self.loops.len() {
            let last = self.last_descendant(i);
            for j in 0..self.loops.len() {
                if i == j {
                    continue;
                }

                let nested = (j > i) && (j <= last);
                const REVERSE_B: bool = false;

                if self.loops[i].contains_non_crossing_boundary(&self.loops[j], REVERSE_B) != nested
                {
                    let nested_str = if !nested { "not " } else { "" };
                    return Err(format!(
                        "invalid nesting: loop {} should {}contain loop {}",
                        i, nested_str, j
                    ));
                }
            }
        }

        Ok(())
    }

    /// intersects_cell reports whether the polygon intersects the given cell.
    pub fn intersects_cell(&self, cell: &Cell) -> bool {
        let mut it = self.index.iterator();
        let relation = it.locate_cell_id(cell.id);

        // If cell does not overlap any index cell, there is no intersection.
        if relation == Disjoint {
            return false;
        }

        // If cell is subdivided into one or more index cells, there is an
        // intersection to within the S2ShapeIndex error bound (see Contains).
        if relation == Subdivided {
            return true;
        }

        // If cell is an index cell, there is an intersection because index cells
        // are created only if they have at least one edge or they are entirely
        // contained by the loop.
        if it.cell_id() == cell.id {
            return true;
        }

        // Otherwise check if any edges intersect cell.
        if self.boundary_approx_intersects(&mut it, cell) {
            return true;
        }

        // Otherwise check if the loop contains the center of cell.
        self.iterator_contains_point(&it, &cell.center())
    }

    /// boundary_approx_intersects reports whether the loop's boundary intersects cell.
    /// It may also return true when the loop boundary does not intersect cell but
    /// some edge comes within the worst-case error tolerance.
    ///
    /// This requires that it.Locate(cell) returned Indexed.
    fn boundary_approx_intersects(&self, it: &mut ShapeIndexIterator, cell: &Cell) -> bool {
        let a_clipped = it
            .index_cell()
            .expect("Why does an indexed cell not exist? Is the ShapeIndex not initialized")
            .find_by_shape_id(0)
            .expect("Somehow there is no shape in this IndexCell!");

        // If there are no edges, there is no intersection.
        if a_clipped.edges.is_empty() {
            return false;
        }

        // We can save some work if cell is the index cell itself.
        if it.cell_id() == cell.id {
            return true;
        }

        // Otherwise check whether any of the edges intersect cell.
        let max_error = FACE_CLIP_ERROR_UV_COORD + INTERSECT_RECT_ERROR_UV_DIST;
        let bound = cell.bound_uv().expanded_by_margin(max_error);

        for e in &a_clipped.edges {
            let edge = self
                .index
                .shape(0)
                .expect("no shape in this index cell at id: 0")
                .edge(*e as i64);
            let uv_coor = clip_to_padded_face(&edge.v0, &edge.v1, cell.face(), max_error);
            if let Some((v0, v1)) = uv_coor
                && edge_intersects_rect(v0, v1, &bound)
            {
                return true;
            }
        }

        false
    }

    /// iterator_contains_point reports whether the iterator that is positioned at the
    /// ShapeIndexCell that may contain p, contains the point p.
    fn iterator_contains_point(&self, it: &ShapeIndexIterator, point: &Point) -> bool {
        // Test containment by drawing a line segment from the cell center to the
        // given point and counting edge crossings.
        let a_clipped = it
            .index_cell()
            .expect("Why does an indexed cell not exist? Is the ShapeIndex not initialized")
            .find_by_shape_id(0)
            .expect("Somehow there is no shape in this IndexCell!");
        let mut inside = a_clipped.contains_center;

        if a_clipped.edges.is_empty() {
            return inside;
        }

        // This block requires ShapeIndex.
        let mut crosser = EdgeCrosser::new(&it.center(), point);
        let shape = self.index.shape(0);

        if let Some(shape) = shape {
            for e in a_clipped.edges.iter() {
                let edge = shape.edge(*e as i64);
                inside = inside != crosser.edge_or_vertex_crossing(&edge.v0, &edge.v1);
            }
        }

        inside
    }

    /// contains_point reports whether the polygon contains the point.
    pub fn contains_point(&self, point: &Point) -> bool {
        // NOTE: A bounds check slows down this function by about 50%. It is
        // worthwhile only when it might allow us to delay building the index.
        if !self.index.is_fresh() && !self.bound.contains_point(point) {
            return false;
        }

        // For small polygons, and during initial construction, it is faster to just
        // check all the crossings.
        const MAX_BRUTE_FORCE_VERTICES: usize = 32;
        if self.num_vertices < MAX_BRUTE_FORCE_VERTICES || self.index.is_empty() {
            let mut inside = false;
            for l in &self.loops {
                // use loops brute force to avoid building the index on each loop.
                inside = inside != l.brute_force_contains_point(point);
            }
            return inside;
        }

        // Otherwise we look up the ShapeIndex cell containing this point.
        // TODO: Implement ContainsPointQuery
        // return NewContainsPointQuery(p.index, VertexModelSemiOpen).Contains(point)
        false // Placeholder
    }
}

impl Region for Polygon {
    // Implement the Region trait

    /// cap_bound returns a bounding spherical cap.
    fn cap_bound(&self) -> Cap {
        self.bound.cap_bound()
    }

    /// rect_bound returns a bounding latitude-longitude rectangle.
    fn rect_bound(&self) -> Rect {
        self.bound.clone()
    }

    /// contains_cell reports whether the polygon contains the given cell.
    fn contains_cell(&self, cell: &Cell) -> bool {
        let mut it = self.index.iterator();
        let relation = it.locate_cell_id(cell.id);

        // If "cell" is disjoint from all index cells, it is not contained.
        // Similarly, if "cell" is subdivided into one or more index cells then it
        // is not contained, since index cells are subdivided only if they (nearly)
        // intersect a sufficient number of edges. (But note that if "cell" itself
        // is an index cell then it may be contained, since it could be a cell with
        // no edges in the loop interior.)
        if relation != Indexed {
            return false;
        }

        // Otherwise check if any edges intersect "cell".
        if self.boundary_approx_intersects(&mut it, cell) {
            return false;
        }

        // Otherwise check if the loop contains the center of "cell".
        self.iterator_contains_point(&it, &cell.center())
    }
}

/// Defines a total ordering on Loops that does not depend on the cyclic
/// order of loop vertices. This function is used to choose which loop to
/// invert in the case where several loops have exactly the same area.
fn compare_loops(a: &Loop, b: &Loop) -> i32 {
    let na = a.num_vertices();
    let nb = b.num_vertices();
    if na != nb {
        return (na as i32) - (nb as i32);
    }

    let (mut ai, a_dir) = a.canonical_first_vertex();
    let (mut bi, b_dir) = b.canonical_first_vertex();

    if a_dir != b_dir {
        let a_dir_as_int = a_dir as i32;
        let b_dir_as_int = b_dir as i32;

        return a_dir_as_int - b_dir_as_int;
    }

    for _n in (0..(a.num_vertices() as i32)).rev() {
        let a_temp = a.vertex(ai + a_dir).0;
        let b_temp = b.vertex(bi + b_dir).0;

        let cmp = a_temp.cmp(&b_temp);
        if cmp != Ordering::Equal {
            match cmp {
                Ordering::Greater => return 1,
                Ordering::Less => return -1,
                _ => {
                    panic!("HOW IS THIS OPSSILBINESRTIENTSROIE! Equal!?!?")
                }
            }
        }

        ai += a_dir;
        bi += b_dir;
    }

    0
}

// Implementation of Shape trait
impl Shape for Polygon {
    /// num_edges returns the number of edges in this shape.
    fn num_edges(&self) -> i64 {
        self.num_edges as i64
    }

    /// edge returns endpoints for the given edge index.
    fn edge(&self, e: i64) -> Edge {
        let e = e as usize;
        let mut i: usize = 0;

        if !self.cumulative_edges.is_empty() {
            for (idx, &_cum_edge) in self.cumulative_edges.iter().enumerate() {
                if idx + 1 >= self.cumulative_edges.len() || e < self.cumulative_edges[idx + 1] {
                    i = idx;
                    break;
                }
            }
        } else {
            // When the number of loops is small, use linear search. Most often
            // there is exactly one loop and the code below executes zero times.
            for (idx, loop_ref) in self.loops.iter().enumerate() {
                if e < loop_ref.num_vertices() {
                    i = idx;
                    break;
                }
            }
        }

        Edge {
            v0: self.loops[i]
                .oriented_vertex(e - self.cumulative_edges.get(i).cloned().unwrap_or(0)),
            v1: self.loops[i]
                .oriented_vertex(e - self.cumulative_edges.get(i).cloned().unwrap_or(0) + 1),
        }
    }

    /// reference_point returns the reference point for this polygon.
    fn reference_point(&self) -> ReferencePoint {
        let mut contains_origin = false;
        for l in &self.loops {
            contains_origin = contains_origin != l.contains_origin();
        }
        ReferencePoint::origin(contains_origin)
    }

    /// num_chains reports the number of contiguous edge chains in the Polygon.
    fn num_chains(&self) -> i64 {
        self.num_loops() as i64
    }

    /// chain returns the i-th edge Chain (loop) in the Shape.
    fn chain(&self, chain_id: i64) -> Chain {
        let chain_id = chain_id as usize;
        if !self.cumulative_edges.is_empty() {
            Chain {
                start: self.cumulative_edges[chain_id] as i64,
                length: self.loops[chain_id].num_vertices() as i64,
            }
        } else {
            let mut e = 0;
            for j in 0..chain_id {
                e += self.loops[j].num_vertices();
            }

            // Polygon represents a full loop as a loop with one vertex, while
            // Shape represents a full loop as a chain with no vertices.
            let num_vertices = self.loops[chain_id].num_vertices();
            if num_vertices != 1 {
                Chain {
                    start: e as i64,
                    length: num_vertices as i64,
                }
            } else {
                Chain {
                    start: e as i64,
                    length: 0,
                }
            }
        }
    }

    /// chain_edge returns the j-th edge of the i-th edge Chain (loop).
    fn chain_edge(&self, i: i64, j: i64) -> Edge {
        let i = i as usize;
        let j = j as usize;
        Edge {
            v0: self.loops[i].oriented_vertex(j),
            v1: self.loops[i].oriented_vertex(j + 1),
        }
    }

    /// chain_position returns a pair (i, j) such that edge_id is the j-th edge
    /// of the i-th edge Chain.
    fn chain_position(&self, edge_id: i64) -> ChainPosition {
        let edge_id = edge_id as usize;
        let mut i: usize = 0;

        if !self.cumulative_edges.is_empty() {
            for (idx, &_cum_edge) in self.cumulative_edges.iter().enumerate() {
                if idx + 1 >= self.cumulative_edges.len()
                    || edge_id < self.cumulative_edges[idx + 1]
                {
                    i = idx;
                    break;
                }
            }
        } else {
            // When the number of loops is small, use linear search. Most often
            // there is exactly one loop and the code below executes zero times.
            for (idx, loop_ref) in self.loops.iter().enumerate() {
                if edge_id < loop_ref.num_vertices() {
                    i = idx;
                    break;
                }
            }
        }

        // TODO: unify this and Edge since they are mostly identical.
        ChainPosition {
            chain_id: i as i64,
            offset: (edge_id - self.cumulative_edges.get(i).cloned().unwrap_or(0)) as i64,
        }
    }

    /// dimension returns the dimension of the geometry represented by this Polygon.
    fn dimension(&self) -> i64 {
        2
    }
}

// Additional geometric methods
impl Polygon {
    /// area returns the area of the polygon interior, i.e. the region on the left side
    /// of an odd number of loops. The return value is between 0 and 4*Pi.
    pub fn area(&self) -> f64 {
        let mut area = 0.0;
        for loop_ref in &self.loops {
            area += loop_ref.sign() as f64 * loop_ref.area();
        }
        area
    }

    /// centroid returns the true centroid of the polygon multiplied by the area of
    /// the polygon. The result is not unit length, so you may want to normalize it.
    /// Also note that in general, the centroid may not be contained by the polygon.
    ///
    /// We prescale by the polygon area for two reasons: (1) it is cheaper to
    /// compute this way, and (2) it makes it easier to compute the centroid of
    /// more complicated shapes (by splitting them into disjoint regions and
    /// adding their centroids).
    pub fn centroid(&self) -> Point {
        let mut u = R3Vector {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        for loop_ref in &self.loops {
            let v = loop_ref.centroid().0;
            if loop_ref.sign() < 0 {
                u = u - v;
            } else {
                u = u + v;
            }
        }
        Point(u)
    }

    /// contains reports whether this polygon contains the other polygon.
    /// Specifically, it reports whether all the points in the other polygon
    /// are also in this polygon.
    pub fn contains(&self, o: &Polygon) -> bool {
        // If both polygons have one loop, use the more efficient Loop method.
        // Note that Loop's Contains does its own bounding rectangle check.
        if self.loops.len() == 1 && o.loops.len() == 1 {
            return self.loops[0].contains(&o.loops[0]);
        }

        // Otherwise if neither polygon has holes, we can still use the more
        // efficient Loop's Contains method (rather than compareBoundary),
        // but it's worthwhile to do our own bounds check first.
        if !self.subregion_bound.contains(&o.bound) {
            // Even though Bound(A) does not contain Bound(B), it is still possible
            // that A contains B. This can only happen when union of the two bounds
            // spans all longitudes. For example, suppose that B consists of two
            // shells with a longitude gap between them, while A consists of one shell
            // that surrounds both shells of B but goes the other way around the
            // sphere (so that it does not intersect the longitude gap).
            if !self.bound.lng.union(&o.bound.lng).is_full() {
                return false;
            }
        }

        if !self.has_holes && !o.has_holes {
            for l in &o.loops {
                if !self.any_loop_contains(l) {
                    return false;
                }
            }
            return true;
        }

        // Polygon A contains B iff B does not intersect the complement of A. From
        // the intersection algorithm below, this means that the complement of A
        // must exclude the entire boundary of B, and B must exclude all shell
        // boundaries of the complement of A. (It can be shown that B must then
        // exclude the entire boundary of the complement of A.) The first call
        // below returns false if the boundaries cross, therefore the second call
        // does not need to check for any crossing edges (which makes it cheaper).
        self.contains_boundary(o) && o.excludes_non_crossing_complement_shells(self)
    }

    /// contains_boundary reports whether this polygon contains the entire boundary of o.
    fn contains_boundary(&self, o: &Polygon) -> bool {
        for l in &o.loops {
            if self.compare_boundary(l) <= 0 {
                return false;
            }
        }
        true
    }

    /// compare_boundary returns +1 if this polygon contains the boundary of o, -1 if A
    /// excludes the boundary of o, and 0 if the boundaries of A and o cross.
    fn compare_boundary(&self, o: &Loop) -> BoundaryCondition {
        let mut result = -1;
        for i in 0..self.loops.len() {
            if result == 0 {
                break;
            }
            // If o crosses any loop of A, the result is 0. Otherwise the result
            // changes sign each time o is contained by a loop of A.

            let boundary_condition_as_int = -self.loops[i].compare_boundary(o) as i32;

            result *= boundary_condition_as_int;
        }
        // TODO: Make this Result return an error.
        result.try_into().unwrap()
    }

    /// excludes_boundary reports whether this polygon excludes the entire boundary of o.
    fn excludes_boundary(&self, o: &Polygon) -> bool {
        for l in &o.loops {
            if self.compare_boundary(l) >= 0 {
                return false;
            }
        }
        true
    }

    /// excludes_non_crossing_shells reports whether given two polygons A and B such that the
    /// boundary of A does not cross any loop of B, if A excludes all shell boundaries of B.
    fn excludes_non_crossing_shells(&self, o: &Polygon) -> bool {
        for l in &o.loops {
            if l.is_hole() {
                continue;
            }
            if self.contains_non_crossing_boundary(l, false) {
                return false;
            }
        }
        true
    }

    /// excludes_non_crossing_complement_shells reports whether given two polygons A and B
    /// such that the boundary of A does not cross any loop of B, if A excludes all
    /// shell boundaries of the complement of B.
    fn excludes_non_crossing_complement_shells(&self, o: &Polygon) -> bool {
        // Special case to handle the complement of the empty or full polygons.
        if o.is_empty() {
            return !self.is_full();
        }
        if o.is_full() {
            return true;
        }

        // Otherwise the complement of B may be obtained by inverting loop(0) and
        // then swapping the shell/hole status of all other loops. This implies
        // that the shells of the complement consist of loop 0 plus all the holes of
        // the original polygon.
        for (j, l) in o.loops.iter().enumerate() {
            if j > 0 && !l.is_hole() {
                continue;
            }

            // The interior of the complement is to the right of loop 0, and to the
            // left of the loops that were originally holes.
            if self.contains_non_crossing_boundary(l, j == 0) {
                return false;
            }
        }
        true
    }

    /// contains_non_crossing_boundary reports whether polygon A contains the boundary of
    /// loop B. Shared edges are handled according to the rule described in loops
    /// contains_non_crossing_boundary.
    fn contains_non_crossing_boundary(&self, o: &Loop, reverse: bool) -> bool {
        let mut inside = false;
        for l in &self.loops {
            let x = l.contains_non_crossing_boundary(o, reverse);
            inside = inside != x;
        }
        inside
    }

    /// intersects reports whether this polygon intersects the other polygon, i.e.
    /// if there is a point that is contained by both polygons.
    pub fn intersects(&self, o: &Polygon) -> bool {
        // If both polygons have one loop, use the more efficient Loop method.
        // Note that Loop Intersects does its own bounding rectangle check.
        if self.loops.len() == 1 && o.loops.len() == 1 {
            return self.loops[0].intersects(&o.loops[0]);
        }

        // Otherwise if neither polygon has holes, we can still use the more
        // efficient Loop.Intersects method. The polygons intersect if and
        // only if some pair of loop regions intersect.
        if !self.bound.intersects(&o.bound) {
            return false;
        }

        if !self.has_holes && !o.has_holes {
            for l in &o.loops {
                if self.any_loop_intersects(l) {
                    return true;
                }
            }
            return false;
        }

        // Polygon A is disjoint from B if A excludes the entire boundary of B and B
        // excludes all shell boundaries of A. (It can be shown that B must then
        // exclude the entire boundary of A.) The first call below returns false if
        // the boundaries cross, therefore the second call does not need to check
        // for crossing edges.
        !self.excludes_boundary(o) || !o.excludes_non_crossing_shells(self)
    }
}
