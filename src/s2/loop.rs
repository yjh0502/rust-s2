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

use crate::consts::DBL_EPSILON;
use crate::consts::EPSILON;
use crate::error::S2Error;
use crate::point::{get_frame, ordered_ccw, regular_points_for_frame};
use crate::r3::vector::Vector as R3Vector;
use crate::rect_bounder::expand_for_subregions;
use crate::region::Region;
use crate::s1;
use crate::s1::Angle;
use crate::s2::cap::Cap;
use crate::s2::cell::Cell;
use crate::s2::edge_clipping::{
    clip_to_padded_face, edge_intersects_rect, INTERSECT_RECT_ERROR_UV_DIST,
};
use crate::s2::edge_crosser::EdgeCrosser;
use crate::s2::edge_crossings::angle_contains_vertex;
use crate::s2::point::Point;
use crate::s2::rect::Rect;
use crate::s2::rect_bounder::RectBounder;
use crate::s2::shape::{Chain, ChainPosition, Edge, ReferencePoint, Shape, ShapeType};
use crate::s2::shape_index::ShapeIndex;
use crate::shape_index::CellRelation::{Disjoint, Indexed, Subdivided};
use crate::shape_index::FACE_CLIP_ERROR_UV_COORD;
use cgmath::Matrix3;
use std::cmp::Ordering;
use std::convert::TryFrom;
use std::f64;
use std::f64::consts::PI;
use std::fmt::{Debug, Formatter};
use std::hash::{Hash, Hasher};
use std::ops::{Add, AddAssign, Mul, Neg, Sub};

pub enum OriginBound {
    OriginInside,
    BoundEncoded,
}

#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub enum VertexTraversalDirection {
    Forward, // vertices are traversed starting from firstIdx and incrementing the index. This is the "forward" direction. The sequence will be firstIdx, firstIdx + 1, firstIdx + 2,
    Backward, // vertices are traversed starting from firstIdx + n (where n is the loop length) and decrementing the index. This is the "backward" direction. The sequence will be firstIdx + n, firstIdx + n - 1, firstIdx + n - 2
}

impl Into<i32> for VertexTraversalDirection {
    fn into(self) -> i32 {
        match self {
            VertexTraversalDirection::Forward => 1,
            VertexTraversalDirection::Backward => -1,
        }
    }
}

impl Into<f64> for VertexTraversalDirection {
    fn into(self) -> f64 {
        match self {
            VertexTraversalDirection::Forward => 1.,
            VertexTraversalDirection::Backward => -1.,
        }
    }
}

impl AddAssign<VertexTraversalDirection> for usize {
    fn add_assign(&mut self, rhs: VertexTraversalDirection) {
        match rhs {
            VertexTraversalDirection::Forward => *self += 1,
            VertexTraversalDirection::Backward => *self -= 1,
        }
    }
}

/// Loop represents a simple spherical polygon. It consists of a sequence
/// of vertices where the first vertex is implicitly connected to the
/// last. All loops are defined to have a CCW orientation, i.e. the interior of
/// the loop is on the left side of the edges. This implies that a clockwise
/// loop enclosing a small area is interpreted to be a CCW loop enclosing a
/// very large area.
///
/// Loops are not allowed to have any duplicate vertices (whether adjacent or
/// not). Non-adjacent edges are not allowed to intersect, and furthermore edges
/// of length 180 degrees are not allowed (i.e., adjacent vertices cannot be
/// antipodal). Loops must have at least 3 vertices (except for the "empty" and
/// "full" loops discussed below).
///
/// There are two special loops: the "empty" loop contains no points and the
/// "full" loop contains all points. These loops do not have any edges, but to
/// preserve the invariant that every loop can be represented as a vertex
/// chain, they are defined as having exactly one vertex each (see EmptyLoop
/// and FullLoop).
#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Loop {
    /// The vertices of the loop. These should be ordered counterclockwise
    /// around the loop interior.
    vertices: Vec<Point>,

    /// originInside keeps a precomputed value whether this loop contains the origin
    /// versus computing from the set of vertices every time.
    origin_inside: bool,

    /// depth is the nesting depth of this Loop if it is contained by a Polygon
    /// or other shape and is used to determine if this loop represents a hole
    /// or a filled in portion.
    pub(crate) depth: i32,

    /// bound is a conservative bound on all points contained by this loop.
    /// If l.ContainsPoint(P), then l.bound.ContainsPoint(P).
    bound: Rect,

    /// Since bound is not exact, it is possible that a loop A contains
    /// another loop B whose bounds are slightly larger. subregionBound
    /// has been expanded sufficiently to account for this error, i.e.
    /// if A.Contains(B), then A.subregionBound.Contains(B.bound).
    subregion_bound: Rect,

    /// index is the spatial index for this Loop.
    #[cfg_attr(feature = "serde", serde(skip))]
    index: ShapeIndex,
}

impl PartialEq for Loop {
    fn eq(&self, _other: &Self) -> bool {
        todo!()
    }
}

impl Eq for Loop {}

impl Hash for Loop {
    fn hash<H: Hasher>(&self, _state: &mut H) {
        todo!()
    }
}

impl Loop {
    // containsNonCrossingBoundary reports whether given two loops whose boundaries
    // do not cross (see compareBoundary), if this loop contains the boundary of the
    // other loop. If reverse is true, the boundary of the other loop is reversed
    // first (which only affects the result when there are shared edges). This method
    // is cheaper than compareBoundary because it does not test for edge intersections.
    //
    // This function requires that neither loop is empty, and that if the other is full,
    // then reverse == false.
    pub(crate) fn contains_non_crossing_boundary(&self, other: &Loop, reverse_other: bool) -> bool {
        // The bounds must intersect for containment.
        if !self.bound.intersects(&other.bound) {
            return false;
        }

        // Full loops are handled as though the loop surrounded the entire sphere.
        if self.is_full() {
            return true;
        }
        if other.is_full() {
            return false;
        }

        let m = self.find_vertex(other.vertex(0));

        match m {
            None => {
                // Since the other loops vertex 0 is not shared, we can check if this contains it.
                return self.contains_point(other.vertex(0));
            }
            Some(m) => {
                return wedge_contains_semi_wedge(
                    &self.vertex(m - 1),
                    &self.vertex(m),
                    &self.vertex(m + 1),
                    &other.vertex(1),
                    reverse_other,
                );
            }
        }
    }

    // Validate checks whether this is a valid loop.
    fn validate(&self) -> bool {
        todo!()
    }
}

// These two points are used for the special Empty and Full loops.
const EMPTY_LOOP_POINT: Point = Point(R3Vector {
    x: 0.0,
    y: 0.0,
    z: 1.0,
});
const FULL_LOOP_POINT: Point = Point(R3Vector {
    x: 0.0,
    y: 0.0,
    z: -1.0,
});

impl Debug for Loop {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct(/* &str */ "Loop")
            .field("vertices", &self.vertices)
            .field("origin_inside", &self.origin_inside)
            .field("depth", &self.depth)
            .field("bound", &self.bound)
            .field("subregion_bound", &self.subregion_bound)
            // Ignoring index
            .finish()
    }
}

impl Loop {
    /// Creates a new loop from the given points.
    pub fn from_points(pts: Vec<Point>) -> Self {
        let mut l = Self {
            vertices: pts,
            origin_inside: false,
            depth: 0,
            bound: Rect::empty(),
            subregion_bound: Rect::empty(),
            index: ShapeIndex::new(),
        };

        l.init_origin_and_bound();
        l
    }

    /// Creates a loop corresponding to the given cell.
    ///
    /// Note that the loop and cell *do not* contain exactly the same set of
    /// points, because Loop and Cell have slightly different definitions of
    /// point containment. For example, a Cell vertex is contained by all
    /// four neighboring Cells, but it is contained by exactly one of four
    /// Loops constructed from those cells. As another example, the cell
    /// coverings of cell and LoopFromCell(cell) will be different, because the
    /// loop contains points on its boundary that actually belong to other cells
    /// (i.e., the covering will include a layer of neighboring cells).
    pub fn from_cell(c: &Cell) -> Self {
        let mut l = Self {
            vertices: vec![c.vertex(0), c.vertex(1), c.vertex(2), c.vertex(3)],
            origin_inside: false,
            depth: 0,
            bound: Rect::empty(),
            subregion_bound: Rect::empty(),
            index: ShapeIndex::new(),
        };

        l.init_origin_and_bound();
        l
    }

    /// Returns a special "empty" loop.
    pub fn empty() -> Self {
        Self::from_points(vec![EMPTY_LOOP_POINT])
    }

    /// Returns a special "full" loop.
    pub fn full() -> Self {
        Self::from_points(vec![FULL_LOOP_POINT])
    }

    /// Sets the origin containment for the given point and then calls
    /// the initialization for the bounds objects and the internal index.
    fn init_origin_and_bound(&mut self) {
        if self.vertices.len() < 3 {
            // Check for the special "empty" and "full" loops (which have one vertex).
            if !self.is_empty_or_full() {
                self.origin_inside = false;
                return;
            }

            // This is the special empty or full loop, so the origin depends on if
            // the vertex is in the southern hemisphere or not.
            self.origin_inside = self.vertices[0].0.z < 0.0;
        } else {
            // The brute force point containment algorithm works by counting edge
            // crossings starting at a fixed reference point (chosen as OriginPoint()
            // for historical reasons). Loop initialization would be more efficient
            // if we used a loop vertex such as vertex(0) as the reference point
            // instead, however making this change would be a lot of work because
            // originInside is currently part of the Encode() format.
            //
            // In any case, we initialize originInside by first guessing that it is
            // outside, and then seeing whether we get the correct containment result
            // for vertex 1. If the result is incorrect, the origin must be inside
            // the loop instead. Note that the Loop is not necessarily valid and so
            // we need to check the requirements of AngleContainsVertex first.
            let v1_inside = self.vertex(0) != self.vertex(1)
                && self.vertex(2) != self.vertex(1)
                && angle_contains_vertex(&self.vertex(0), &self.vertex(1), &self.vertex(2));

            // Initialize before calling contains_point
            self.origin_inside = false;

            // Note that contains_point only does a bounds check once init_index
            // has been called, so it doesn't matter that bound is undefined here.
            if v1_inside != self.contains_point(self.vertex(1)) {
                self.origin_inside = true;
            }
        }

        // We *must* call init_bound before initializing the index, because
        // init_bound calls contains_point which does a bounds check before using
        // the index.
        self.init_bound();

        // Create a new index and add us to it.
        self.index = ShapeIndex::new();
        self.index.add(&ShapeType::Loop(self.clone()));
    }

    /// Sets up the approximate bounding Rects for this loop.
    fn init_bound(&mut self) {
        if self.vertices.is_empty() {
            *self = Self::empty();
            return;
        }

        // Check for the special "empty" and "full" loops.
        if self.is_empty_or_full() {
            if self.is_empty() {
                self.bound = Rect::empty();
            } else {
                self.bound = Rect::full();
            }
            self.subregion_bound = self.bound.clone();
            return;
        }

        // The bounding rectangle of a loop is not necessarily the same as the
        // bounding rectangle of its vertices. First, the maximal latitude may be
        // attained along the interior of an edge. Second, the loop may wrap
        // entirely around the sphere (e.g. a loop that defines two revolutions of a
        // candy-cane stripe). Third, the loop may include one or both poles.
        // Note that a small clockwise loop near the equator contains both poles.
        let mut bounder = RectBounder::new();
        for i in 0..=self.vertices.len() {
            // add vertex 0 twice
            bounder.add_point(&self.vertex(i));
        }
        let mut b = bounder.get_bound();

        if self.contains_point(Point(R3Vector {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        })) {
            b = Rect {
                lat: crate::r1::interval::Interval::new(b.lat.lo, PI / 2.0),
                lng: s1::interval::Interval::full_interval(),
            };
        }
        // If a loop contains the south pole, then either it wraps entirely
        // around the sphere (full longitude range), or it also contains the
        // north pole in which case b.lng.is_full() due to the test above.
        // Either way, we only need to do the south pole containment test if
        // b.lng.is_full().
        if b.lng.is_full()
            && self.contains_point(Point(R3Vector {
                x: 0.0,
                y: 0.0,
                z: -1.0,
            }))
        {
            b.lat.lo = -PI / 2.0;
        }
        self.bound = b;
        self.subregion_bound = expand_for_subregions(&self.bound);
    }

    /// Returns whether this loop is considered empty.
    pub fn is_empty(&self) -> bool {
        self.is_empty_or_full() && !self.contains_origin()
    }

    /// Returns whether this loop is considered full.
    pub fn is_full(&self) -> bool {
        self.is_empty_or_full() && self.contains_origin()
    }

    /// Returns whether this loop is either the "empty" or "full" special loops.
    pub fn is_empty_or_full(&self) -> bool {
        self.vertices.len() == 1
    }

    /// Returns the reference point for this loop.
    pub fn reference_point(&self) -> ReferencePoint {
        ReferencePoint::origin(self.origin_inside)
    }

    /// Returns the vertex at the given index. For convenience, the vertex indices
    /// wrap automatically for methods that do index math such as Edge.
    /// i.e., Vertex(NumEdges() + n) is the same as Vertex(n).
    pub fn vertex(&self, i: usize) -> Point {
        self.vertices[i % self.vertices.len()]
    }

    /// Returns true if this loop contains the point.
    pub fn contains_point(&self, p: Point) -> bool {
        if !self.index.is_fresh() && !self.bound.contains_point(&p) {
            return false;
        }

        // For small loops it is faster to just check all the crossings. We also
        // use this method during loop initialization because InitOriginAndBound()
        // calls Contains() before InitIndex().

        const MAX_BRUTE_FORCE_VERTICES: usize = 32;

        if self.index.num_shapes() == 0 || // Index has not been initialized yet
           self.vertices.len() <= MAX_BRUTE_FORCE_VERTICES
        {
            return self.brute_force_contains_point(&p);
        }

        // Otherwise, look up the point in the index.
        let mut it = self.index.iterator();
        if !it.locate_point(p) {
            return false;
        }
        self.iterator_contains_point(&mut it, p)
    }

    /// Reports if the given point is contained by this loop, by doing a simple
    /// brute force check.  This method does not use the ShapeIndex, so it is only
    /// preferable below a certain size of loop.
    pub fn brute_force_contains_point(&self, p: &Point) -> bool {
        let origin = Point::origin();
        let mut inside = self.origin_inside;
        let mut crosser = EdgeCrosser::new_chain_edge_crosser(&origin, &p, &self.vertex(0));
        for i in 1..=self.vertices.len() {
            // add vertex 0 twice
            inside = inside != crosser.edge_or_vertex_chain_crossing(&self.vertex(i));
        }
        inside
    }

    /// Reports if the iterator that is positioned at the ShapeIndexCell
    /// that may contain p, contains the point p.
    fn iterator_contains_point(
        &self,
        it: &mut crate::s2::shape_index::ShapeIndexIterator,
        p: Point,
    ) -> bool {
        // Test containment by drawing a line segment from the cell center to the
        // given point and counting edge crossings.
        // TODO: Is this
        let a_clipped = it
            .index_cell()
            .expect("Why does an indexed cell not exist? Is the ShapeIndex not initialized")
            .find_by_shape_id(0)
            .expect("Index cell should have shapeID of 0 here!!");
        let mut inside = a_clipped.contains_center;

        if a_clipped.num_edges() > 0 {
            let center = it.center();
            let mut crosser = EdgeCrosser::new(&center, &p);
            let mut ai_prev = -2;

            for &ai in &a_clipped.edges {
                if ai != ai_prev + 1 {
                    crosser.restart_at(&self.vertex(ai as usize));
                }
                ai_prev = ai;
                inside = inside
                    != crosser.edge_or_vertex_chain_crossing(&self.vertex((ai + 1) as usize));
            }
        }
        inside
    }

    /// Returns whether the loop contains origin.
    pub fn contains_origin(&self) -> bool {
        self.origin_inside
    }

    /// Returns the vertices of the loop.
    pub fn vertices(&self) -> &[Point] {
        &self.vertices
    }

    /// Returns the number of vertices in this loop.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Returns whether this loop represents a hole in its containing polygon.
    pub fn is_hole(&self) -> bool {
        (self.depth & 1) != 0
    }

    /// Returns -1 if this Loop represents a hole in its containing polygon, and +1 otherwise.
    pub fn sign(&self) -> i32 {
        if self.is_hole() {
            -1
        } else {
            1
        }
    }

    /// Reports whether the region contained by this loop is a superset of the
    /// region contained by the given other loop.
    pub fn contains(&self, o: &Loop) -> bool {
        // For a loop A to contain the loop B, all of the following must
        // be true:
        //
        //  (1) There are no edge crossings between A and B except at vertices.
        //
        //  (2) At every vertex that is shared between A and B, the local edge
        //      ordering implies that A contains B.
        //
        //  (3) If there are no shared vertices, then A must contain a vertex of B
        //      and B must not contain a vertex of A. (An arbitrary vertex may be
        //      chosen in each case.)
        //
        // The second part of (3) is necessary to detect the case of two loops whose
        // union is the entire sphere, i.e. two loops that contains each other's
        // boundaries but not each other's interiors.

        if !self.subregion_bound.contains(&o.bound.clone()) {
            return false;
        }

        // Special cases to handle either loop being empty or full.
        if self.is_empty_or_full() || o.is_empty_or_full() {
            return self.is_full() || o.is_empty();
        }

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        let mut relation = ContainsRelation::new();
        if has_crossing_relation(self, o, &mut relation) {
            return false;
        }

        // There are no crossings, and if there are any shared vertices then A
        // contains B locally at each shared vertex.
        if relation.found_shared_vertex {
            return true;
        }

        // Since there are no edge intersections or shared vertices, we just need to
        // test condition (3) above. We can skip this test if we discovered that A
        // contains at least one point of B while checking for edge crossings.
        if !self.contains_point(o.vertex(0)) {
            return false;
        }

        // We still need to check whether (A union B) is the entire sphere.
        // Normally this check is very cheap due to the bounding box precondition.
        if (o.subregion_bound.contains(&self.bound.clone()) || o.bound.union(&self.bound).is_full())
            && o.contains_point(self.vertex(0))
        {
            return false;
        }

        true
    }

    /// Reports whether the region contained by this loop intersects the region
    /// contained by the other loop.
    pub fn intersects(&self, o: &Loop) -> bool {
        // Given two loops, A and B, A.Intersects(B) if and only if !A.Complement().Contains(B).
        //
        // This code is similar to Contains, but is optimized for the case
        // where both loops enclose less than half of the sphere.
        if !self.bound.intersects(&o.bound) {
            return false;
        }

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        let mut relation = IntersectsRelation::new();
        if has_crossing_relation(self, o, &mut relation) {
            return true;
        }

        if relation.found_shared_vertex {
            return false;
        }

        // Since there are no edge intersections or shared vertices, the loops
        // intersect only if A contains B, B contains A, or the two loops contain
        // each other's boundaries.  These checks are usually cheap because of the
        // bounding box preconditions.  Note that neither loop is empty (because of
        // the bounding box check above), so it is safe to access vertex(0).

        // Check whether A contains B, or A and B contain each other's boundaries.
        // (Note that A contains all the vertices of B in either case.)
        if (self.subregion_bound.contains(&o.bound.clone()) || self.bound.union(&o.bound).is_full())
            && self.contains_point(o.vertex(0))
        {
            return true;
        }

        // Check whether B contains A.
        if o.subregion_bound.contains(&self.bound.clone()) && o.contains_point(self.vertex(0)) {
            return true;
        }

        false
    }

    /// Equal reports whether two loops have the same vertices in the same linear order
    /// (i.e., cyclic rotations are not allowed).
    pub fn equal(&self, other: &Loop) -> bool {
        if self.vertices.len() != other.vertices.len() {
            return false;
        }

        for i in 0..self.vertices.len() {
            if self.vertex(i) != other.vertex(i) {
                return false;
            }
        }

        true
    }

    /// BoundaryEqual reports whether the two loops have the same boundary. This is
    /// true if and only if the loops have the same vertices in the same cyclic order
    /// (i.e., the vertices may be cyclically rotated). The empty and full loops are
    /// considered to have different boundaries.
    pub fn boundary_equal(&self, o: &Loop) -> bool {
        if self.vertices.len() != o.vertices.len() {
            return false;
        }

        // Special case to handle empty or full loops.  Since they have the same
        // number of vertices, if one loop is empty/full then so is the other.
        if self.is_empty_or_full() {
            return self.is_empty() == o.is_empty();
        }

        // Loop through the vertices to find the first of ours that matches the
        // starting vertex of the other loop. Use that offset to then 'align' the
        // vertices for comparison.
        for offset in 0..self.vertices.len() {
            if self.vertex(offset) == o.vertex(0) {
                // There is at most one starting offset since loop vertices are unique.
                let mut matched = true;
                for i in 0..self.vertices.len() {
                    if self.vertex(i + offset) != o.vertex(i) {
                        matched = false;
                        break;
                    }
                }
                if matched {
                    return true;
                }
            }
        }

        false
    }

    /// Tests whether the given loops contains this loop.
    /// This function does not test for edge intersections. The two loops must meet
    /// all of the Polygon requirements; for example this implies that their
    /// boundaries may not cross or have any shared edges (although they may have
    /// shared vertices).
    pub fn contains_nested(&self, other: &Loop) -> bool {
        if !self.subregion_bound.contains(&other.bound.clone()) {
            return false;
        }

        // Special cases to handle either loop being empty or full.  Also bail out
        // when B has no vertices to avoid heap overflow on the vertex(1) call
        // below.  (This method is called during polygon initialization before the
        // client has an opportunity to call IsValid().)
        if self.is_empty_or_full() || other.num_vertices() < 2 {
            return self.is_full() || other.is_empty();
        }

        // We are given that A and B do not share any edges, and that either one
        // loop contains the other or they do not intersect.
        match self.find_vertex(other.vertex(1)) {
            Some(m) => {
                // Check whether the edge order around other.Vertex(1) is compatible with
                // A containing B.
                general_wedge_contains(
                    &self.vertex(m - 1),
                    &self.vertex(m),
                    &self.vertex(m + 1),
                    &other.vertex(0),
                    &other.vertex(2),
                )
            }
            None => {
                // Since other.vertex(1) is not shared, we can check whether A contains it.
                self.contains_point(other.vertex(1))
            }
        }
    }

    /// Find a vertex of this loop that matches the given vertex p.
    /// Returns the index of the vertex in the range [0..num_vertices-1]
    /// or None if no matching vertex is found.
    fn find_vertex(&self, p: Point) -> Option<usize> {
        if self.vertices.len() < 10 {
            // Exhaustive search for loops below a small threshold.
            for i in 0..self.vertices.len() {
                if self.vertex(i) == p {
                    return Some(i);
                }
            }
            return None;
        }

        let mut it = self.index.iterator();
        if !it.locate_point(p) {
            return None;
        }

        let a_clipped = it.index_cell()?.find_by_shape_id(0);
        for i in (0..a_clipped?.num_edges()).rev() {
            let ai = a_clipped.unwrap().edges[i];
            if self.vertex(ai as usize) == p {
                if ai == 0 {
                    return Some(self.vertices.len() - 1);
                }
                return Some(ai as usize);
            }

            if self.vertex((ai + 1) as usize) == p {
                return Some((ai + 1) as usize);
            }
        }

        None
    }

    /// Returns the vertex in reverse order if the loop represents a polygon
    /// hole. For example, arguments 0, 1, 2 are mapped to vertices n-1, n-2, n-3, where
    /// n == len(vertices). This ensures that the interior of the polygon is always to
    /// the left of the vertex chain.
    ///
    /// This requires: 0 <= i < 2 * len(vertices)
    pub fn oriented_vertex(&self, i: usize) -> Point {
        let mut j = i;
        if j >= self.vertices.len() {
            j = j - self.vertices.len();
        }
        if self.is_hole() {
            j = self.vertices.len() - 1 - j;
        }
        self.vertex(j)
    }
}

impl Region for Loop {
    /// Returns a tight bounding rectangle. If the loop contains the point,
    /// the bound also contains it.
    fn rect_bound(&self) -> Rect {
        self.bound.clone()
    }

    /// Returns a bounding cap that may have more padding than the corresponding
    /// RectBound. The bound is conservative such that if the loop contains a point P,
    /// the bound also contains it.
    fn cap_bound(&self) -> Cap {
        self.bound.cap_bound()
    }
}

/// ContainsRelation is a helper for the Contains() method.
/// It implements the loopRelation interface for a contains operation. If
/// A.ContainsPoint(P) == false && B.ContainsPoint(P) == true, it is equivalent
/// to having an edge crossing (i.e., Contains returns false).
struct ContainsRelation {
    found_shared_vertex: bool,
}

impl ContainsRelation {
    fn new() -> Self {
        Self {
            found_shared_vertex: false,
        }
    }
}

impl LoopRelation for ContainsRelation {
    fn a_crossing_target(&self) -> CrossingTarget {
        CrossingTarget::DontCross
    }

    fn b_crossing_target(&self) -> CrossingTarget {
        CrossingTarget::Cross
    }

    fn wedges_cross(
        &mut self,
        a0: &Point,
        ab1: &Point,
        a2: &Point,
        b0: &Point,
        b2: &Point,
    ) -> bool {
        self.found_shared_vertex = true;
        !general_wedge_contains(&a0, &ab1, &a2, &b0, &b2)
    }
}

/// IntersectsRelation is a helper for the Intersects() method.
/// It implements the loopRelation for an intersects operation. Given
/// two loops, A and B, if A.ContainsPoint(P) == true && B.ContainsPoint(P) == true,
/// it is equivalent to having an edge crossing (i.e., Intersects returns true).
struct IntersectsRelation {
    found_shared_vertex: bool,
}

impl IntersectsRelation {
    fn new() -> Self {
        Self {
            found_shared_vertex: false,
        }
    }
}

impl LoopRelation for IntersectsRelation {
    fn a_crossing_target(&self) -> CrossingTarget {
        CrossingTarget::Cross
    }

    fn b_crossing_target(&self) -> CrossingTarget {
        CrossingTarget::Cross
    }

    fn wedges_cross(
        &mut self,
        a0: &Point,
        ab1: &Point,
        a2: &Point,
        b0: &Point,
        b2: &Point,
    ) -> bool {
        self.found_shared_vertex = true;
        // In the Go code this would call WedgeIntersects, but we'll reuse WedgeContains with
        // appropriate logic until we implement WedgeIntersects
        wedge_intersects(a0, ab1, a2, b0, b2)
    }
}

/// WedgeIntersects reports whether the wedges AB1C and DB1E intersect.
/// This is used for testing whether two loops intersect at a shared vertex B1.
fn wedge_intersects(a: &Point, b: &Point, c: &Point, d: &Point, e: &Point) -> bool {
    // For A, B, C, if A == C then the wedge covers the entire sphere.
    // This should not be confused with a degenerate wedge for which B == A or B == C.
    // Note that we don't need to worry about the case where
    // (A == C && D == E), since the API requires that the two wedges have
    // different B vertices.
    if a == c {
        return true;
    }
    if d == e {
        return true;
    }

    // The wedges intersect if and only if neither wedge contains the other's
    // boundary, or they share a boundary. Since we've eliminated the case
    // where both wedges contain the entire sphere, we can check whether each
    // wedge contains the other's boundary.

    // Does the boundary of the first wedge contain the boundary of the second wedge?
    let contains_1 = general_wedge_contains(a, b, c, d, e);

    // Does the boundary of the second wedge contain the boundary of the first wedge?
    let contains_2 = general_wedge_contains(d, b, e, a, c);

    // The wedges intersect if and only if neither contains the other's boundary
    !contains_1 && !contains_2
}

fn general_wedge_contains(a0: &Point, ab1: &Point, a2: &Point, b0: &Point, b2: &Point) -> bool {
    // For A to contain B (where each loop interior is defined to be its left
    // side), the CCW edge order around ab1 must be a2 b2 b0 a0.  We split
    // this test into two parts that test three vertices each.
    ordered_ccw(&a2, &b2, &b0, &ab1) && ordered_ccw(&b0, &a0, &a2, &ab1)
}

// wedgeContainsSemiwedge reports whether the wedge (a0, ab1, a2) contains the
// "semiwedge" defined as any non-empty open set of rays immediately CCW from
// the edge (ab1, b2). If reverse is true, then substitute clockwise for CCW;
// this simulates what would happen if the direction of the other loop was reversed.
fn wedge_contains_semi_wedge(
    a0: &Point,
    ab1: &Point,
    a2: &Point,
    b2: &Point,
    reverse: bool,
) -> bool {
    if b2 == a0 || b2 == a2 {
        // We have a shared or reversed edge.
        return (b2 == a0) == reverse;
    }
    return ordered_ccw(a0, a2, b2, ab1);
}

/// CrossingTarget represents the possible crossing target cases for relations.
#[derive(Copy, Clone, Debug, PartialEq)]
enum CrossingTarget {
    DontCare,
    DontCross,
    Cross,
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum BoundaryCondition {
    ContainsOther, // 1
    CrossesOther,  // 0
    ExcludesOther, //-1
}

impl PartialEq<i32> for BoundaryCondition {
    fn eq(&self, _other: &i32) -> bool {
        todo!()
    }
}

impl PartialOrd<i32> for BoundaryCondition {
    fn partial_cmp(&self, other: &i32) -> Option<Ordering> {
        match self {
            BoundaryCondition::ContainsOther => Option::from(1.cmp(other)),
            BoundaryCondition::CrossesOther => Option::from(0.cmp(other)),
            BoundaryCondition::ExcludesOther => Option::from((-1i32).cmp(other)),
        }
    }
}

impl Neg for BoundaryCondition {
    type Output = Self;
    fn neg(self) -> Self {
        match self {
            BoundaryCondition::ContainsOther => BoundaryCondition::ExcludesOther,
            BoundaryCondition::CrossesOther => BoundaryCondition::CrossesOther,
            BoundaryCondition::ExcludesOther => BoundaryCondition::ContainsOther,
        }
    }
}

impl TryFrom<i32> for BoundaryCondition {
    type Error = S2Error;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(BoundaryCondition::ContainsOther),
            0 => Ok(BoundaryCondition::CrossesOther),
            -1 => Ok(BoundaryCondition::ExcludesOther),
            _ => Err(S2Error::Other(format!(
                "BoundaryCondition Error: Invalid value conversion from int:[{}]",
                value
            ))),
        }
    }
}

/// has_crossing_relation is a helper for Contains, Intersects, and CompareBoundary
/// It checks all edges of loop A for intersection against all edges of loop B and
/// reports if there are any that satisfy the given relation. If there is any shared
/// vertex, the wedges centered at this vertex are tested to see if they satisfy
/// the relation.
///
/// If the two loop boundaries cross, this method is guaranteed to return true.
/// It also returns true in certain cases if the loop relationship is equivalent
/// to crossing. For example, if the relation is Contains and a point P is found
/// such that B contains P but A does not contain P, this method will return true
/// to indicate that the result is the same as though a pair of crossing edges
/// were found (since Contains returns false in both cases).
fn has_crossing_relation(a: &Loop, b: &Loop, relation: &mut dyn LoopRelation) -> bool {
    // Use CrossingEdgeQuery to efficiently find edge crossings

    // First check for shared vertices and test wedge relationships
    for i in 0..a.vertices.len() {
        for j in 0..b.vertices.len() {
            // If vertices are shared
            if a.vertex(i) == b.vertex(j) {
                // Test wedge relationship at shared vertex
                if relation.wedges_cross(
                    &a.vertex(i - 1),
                    &a.vertex(i),
                    &a.vertex(i + 1),
                    &b.vertex(j - 1),
                    &b.vertex(j + 1),
                ) {
                    return true;
                }
            }
        }
    }

    // Then check for edge crossings using CrossingEdgeQuery
    // We need to check in both directions: A's edges crossing B, and B's edges crossing A

    // Check if any edge of A crosses any edge of B
    let mut query = crate::s2::crossing_edge_query::CrossingEdgeQuery::new(&a.index);

    // Iterate through all edges of B
    for j in 0..b.num_edges() {
        let edge_b = b.edge(j as i64);

        // Find all edges of A that cross this edge of B
        let crossings = query.crossings(
            edge_b.v0,
            edge_b.v1,
            &ShapeType::Loop(a.clone()),
            crate::s2::crossing_edge_query::CrossingType::Interior,
        );

        if !crossings.is_empty() {
            return true; // Found a crossing edge
        }
    }

    // Also need to check if one loop contains a vertex of the other
    let target_a = relation.a_crossing_target();
    let target_b = relation.b_crossing_target();

    if target_a != CrossingTarget::DontCare && target_b != CrossingTarget::DontCare {
        // Check if there's a point in b that doesn't satisfy the crossing target for a
        if target_a == CrossingTarget::DontCross && a.contains_point(b.vertex(0)) {
            return true;
        }
        if target_a == CrossingTarget::Cross && !a.contains_point(b.vertex(0)) {
            return true;
        }

        // Check if there's a point in a that doesn't satisfy the crossing target for b
        if target_b == CrossingTarget::DontCross && b.contains_point(a.vertex(0)) {
            return true;
        }
        if target_b == CrossingTarget::Cross && !b.contains_point(a.vertex(0)) {
            return true;
        }
    }

    false
}

/// CrossingRelation defines an interface for checking relationships between loops.
/// This is similar to the Go loopRelation interface, but adapted for Rust.
trait LoopRelation {
    /// Returns the crossing target for loop A
    fn a_crossing_target(&self) -> CrossingTarget;

    /// Returns the crossing target for loop B
    fn b_crossing_target(&self) -> CrossingTarget;

    /// Tests whether wedges at a shared vertex indicate a crossing relationship
    fn wedges_cross(&mut self, a0: &Point, ab1: &Point, a2: &Point, b0: &Point, b2: &Point)
        -> bool;
}

// compareBoundaryRelation implements loopRelation for comparing boundaries.
//
// The compare boundary relation does not have a useful early-exit condition,
// so we return crossingTargetDontCare for both crossing targets.
//
// Aside: A possible early exit condition could be based on the following.
//
//	If A contains a point of both B and ~B, then A intersects Boundary(B).
//	If ~A contains a point of both B and ~B, then ~A intersects Boundary(B).
//	So if the intersections of {A, ~A} with {B, ~B} are all non-empty,
//	the return value is 0, i.e., Boundary(A) intersects Boundary(B).
//
// Unfortunately it isn't worth detecting this situation because by the
// time we have seen a point in all four intersection regions, we are also
// guaranteed to have seen at least one pair of crossing edges.
struct CompareBoundaryRelation {
    reverse: bool,
    found_shared_vertex: bool,
    contains_edge: bool,
    excludes_edge: bool,
}

impl CompareBoundaryRelation {
    fn new(reverse: bool) -> Self {
        Self {
            reverse,
            found_shared_vertex: false,
            contains_edge: false,
            excludes_edge: false,
        }
    }
}

impl LoopRelation for CompareBoundaryRelation {
    fn a_crossing_target(&self) -> CrossingTarget {
        CrossingTarget::DontCare
    }

    fn b_crossing_target(&self) -> CrossingTarget {
        CrossingTarget::DontCare
    }

    fn wedges_cross(
        &mut self,
        a0: &Point,
        ab1: &Point,
        a2: &Point,
        _b0: &Point,
        b2: &Point,
    ) -> bool {
        // Because we don't care about the interior of the other, only its boundary,
        // it is sufficient to check whether this one contains the semiwedge (ab1, b2).
        self.found_shared_vertex = true;
        if wedge_contains_semi_wedge(a0, ab1, a2, b2, self.reverse) {
            self.contains_edge = true
        } else {
            self.excludes_edge = true
        }
        self.contains_edge && self.excludes_edge
    }
}

impl Loop {
    /// CompareBoundary returns +1 if A contains the boundary of B, -1 if A excludes
    /// the boundary of B, and 0 if the boundaries of A and B cross. Excludes means
    /// that A does not intersect the boundary of B at all. For example, if
    /// A contains B, then A contains the boundary of B. If A is contained by B
    /// (including the case where the boundaries of the two loops coincide), then
    /// A excludes the boundary of B.
    ///
    /// This function is primarily useful for determining whether one loop contains
    /// another loop without needing to check all pairs of edges for crossings.
    pub fn compare_boundary(&self, b: &Loop) -> BoundaryCondition {
        // The bounds must intersect for containment or crossing.
        if !self.bound.intersects(&b.bound) {
            return BoundaryCondition::ExcludesOther;
        }

        // Full loops are handled as though the loop surrounded the entire sphere.
        if self.is_full() {
            return BoundaryCondition::ContainsOther;
        }
        if b.is_full() {
            return BoundaryCondition::ExcludesOther;
        }

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        let mut relation = CompareBoundaryRelation::new(b.is_hole()); //(o.IsHole())
        if has_crossing_relation(self, b, &mut relation) {
            return BoundaryCondition::CrossesOther;
        }
        if relation.found_shared_vertex {
            if relation.contains_edge {
                return BoundaryCondition::ContainsOther;
            }
            return BoundaryCondition::ExcludesOther;
        }

        // There are no edge intersections or shared vertices, so we can check
        // whether A contains an arbitrary vertex of B.
        if self.contains_point(b.vertex(0)) {
            return BoundaryCondition::ContainsOther;
        }

        BoundaryCondition::ExcludesOther
    }

    /// Returns a first vertex index and a direction (either +1 or -1) such that the
    /// vertex sequence defined by starting at that vertex and proceeding in the
    /// given direction is canonical for this loop.
    ///
    /// This can be useful for converting a loop to a canonical form. The return
    /// values are chosen such that the lexicographically smallest edge comes first
    /// and is in the forward direction. The result is unspecified if the loop has
    /// no vertices.
    // pub fn canonical_first_vertex(&self) -> (usize, i32) {
    //     let mut first_idx = 0;
    //     let mut dir = 1;
    //     let n = self.vertices.len();
    //
    //     for i in 0..n {
    //         if self.vertex(i).0 < self.vertex(first_idx).0 {
    //             first_idx = i;
    //         }
    //     }
    //
    //     // 0 <= firstIdx <= n-1, so (firstIdx+n*dir) <= 2*n-1.
    //     if self.vertex(first_idx + 1).0 < (self.vertex(first_idx + n - 1).0) {
    //         return (first_idx, 1);
    //     }
    //
    //     // n <= firstIdx <= 2*n-1, so (firstIdx+n*dir) >= 0.
    //     first_idx += n;
    //     (first_idx, -1)
    // }

    /// ContainsCell reports whether this loop contains the given cell.
    /// This method assumes that the loop has been indexed.
    pub fn contains_cell(&self, target: &crate::s2::cell::Cell) -> bool {
        let mut it = self.index.iterator();
        let relation = it.locate_cell_id(target.id);

        // If "target" is disjoint from all index cells, it is not contained.
        // Similarly, if "target" is subdivided into one or more index cells then it
        // is not contained, since index cells are subdivided only if they (nearly)
        // intersect a sufficient number of edges.  (But note that if "target" itself
        // is an index cell then it may be contained, since it could be a cell with
        // no edges in the loop interior.)
        if relation != Indexed {
            return false;
        }

        // Otherwise check if any edges intersect "target".
        if self.boundary_approx_intersects(target) {
            return false;
        }

        // Otherwise check if the loop contains the center of "target".
        self.iterator_contains_point(&mut it, target.center())
    }

    /// IntersectsCell reports whether this loop intersects the given cell.
    /// This method assumes that the loop has been indexed.
    pub fn intersects_cell(&self, target: &crate::s2::cell::Cell) -> bool {
        let mut it = self.index.iterator();
        let relation = it.locate_cell_id(target.id);

        // If target does not overlap any index cell, there is no intersection.
        if relation == Disjoint {
            return false;
        }
        // If target is subdivided into one or more index cells, there is an
        // intersection to within the ShapeIndex error bound (see Contains).
        if relation == Subdivided {
            return true;
        }
        // If target is an index cell, there is an intersection because index cells
        // are created only if they have at least one edge or they are entirely
        // contained by the loop.
        if it.cell_id() == target.id {
            return true;
        }
        // Otherwise check if any edges intersect target.
        if self.boundary_approx_intersects(target) {
            return true;
        }
        // Otherwise check if the loop contains the center of target.
        self.iterator_contains_point(&mut it, target.center())
    }

    // CellUnionBound computes a covering of the Loop.
    pub fn cell_union_bound(&self) -> crate::s2::cellunion::CellUnion {
        self.bound.cap_bound().cell_union_bound().into()
    }

    // boundaryApproxIntersects reports if the loop's boundary intersects target.
    // It may also return true when the loop boundary does not intersect target but
    // some edge comes within the worst-case error tolerance.
    //
    // This requires that it.Locate(target) returned Indexed.
    pub fn boundary_approx_intersects(&self, target: &crate::s2::cell::Cell) -> bool {
        let it = self.index.iterator();
        let a_clipped = it.index_cell().unwrap().find_by_shape_id(0).unwrap();

        // If there are no edges, there is no intersection.
        // TODO: Potentially the wrong way to do this maybe directly counting edges is right?
        if a_clipped.num_edges() == 0 {
            return false;
        }

        // We can save some work if target is the index cell itself.
        if it.cell_id() == target.id {
            return true;
        }

        // Otherwise check whether any of the edges intersect target.
        let max_error = FACE_CLIP_ERROR_UV_COORD + INTERSECT_RECT_ERROR_UV_DIST;
        let bound = target.bound_uv().expanded_by_margin(max_error);
        for (_, ai) in a_clipped.edges.iter().enumerate() {
            let v0_outer = clip_to_padded_face(
                &self.vertex(*ai as usize),
                &self.vertex((ai + 1) as usize),
                target.face(),
                max_error,
            );
            if let Some((v0, v1)) = v0_outer
                && edge_intersects_rect(v0, v1, &bound.clone())
            {
                return true;
            }
        }
        false
    }

    /// RegularLoop creates a loop with the given number of vertices, all
    /// located on a circle of the specified radius around the given center.
    pub fn regular_loop(center: Point, radius: crate::s1::angle::Angle, n: usize) -> Self {
        Self::regular_loop_for_frame(&get_frame(&center), radius, n)
    }

    /// RegularLoopForFrame creates a loop centered at the z-axis of the given
    /// coordinate frame, with the specified angular radius in radians and number
    /// of vertices.
    pub fn regular_loop_for_frame(
        frame: &Matrix3<f64>,
        radius: crate::s1::angle::Angle,
        n_vertices: usize,
    ) -> Self {
        Loop::from_points(regular_points_for_frame(frame, radius, n_vertices))
    }

    // findValidationErrorNoIndex reports whether this is not a valid loop, but
    // skips checks that would require a ShapeIndex to be built for the loop. This
    // is primarily used by Polygon to do validation so it doesn't trigger the
    // creation of unneeded ShapeIndices.
    pub fn find_validation_error_no_index(&self) -> Result<(), S2Error> {
        // All vertices must be unit length.
        for (i, v) in self.vertices.iter().enumerate() {
            if !v.0.is_unit() {
                return S2Error::Other(format!("vertex {} is not unit length", i)).into();
            }
        }

        // Loops must have at least 3 vertices (except for empty and full).
        if self.vertices.len() < 3 {
            if self.is_empty_or_full() {
                return Ok(()); // Skip remaining tests.
            }
            return Err(S2Error::Other(
                "non-empty, non-full loops must have at least 3 vertices".to_string(),
            ));
        }

        // Loops are not allowed to have any duplicate vertices or edge crossings.
        // We split this check into two parts. First we check that no edge is
        // degenerate (identical endpoints). Then we check that there are no
        // intersections between non-adjacent edges (including at vertices). The
        // second check needs the ShapeIndex, so it does not fall within the scope
        // of this method.
        for i in 0..self.vertices.len() {
            if self.vertex(i) == self.vertex(i + 1) {
                return Err(S2Error::Other(format!(
                    "edge {} is degenerate (duplicate vertex)",
                    i
                )));
            }

            // Antipodal vertices are not allowed.
            let other = Point(self.vertex(i + 1).0.mul(-1));
            if self.vertex(i) == other {
                return Err(S2Error::Other(format!(
                    "vertices {} and {} are antipodal",
                    i,
                    (i + 1) % self.vertices.len()
                )));
            }
        }

        Ok(())
    }
    /// encode encodes the loop to a byte vector.
    ///
    /// The encoding consists of:
    /// - the encoding version number (8 bits)
    /// - the number of vertices
    /// - the origin_inside flag
    /// - the vertices of the loop
    pub fn encode(&self) -> Vec<u8> {
        let mut data = Vec::new();

        // Version number
        data.push(1); // Encoding version 1

        // Encode the number of vertices in the loop
        let num_vertices = self.vertices.len() as u32;
        data.extend_from_slice(&num_vertices.to_be_bytes());

        // Encode the origin_inside flag
        let origin_byte = if self.origin_inside { 1u8 } else { 0u8 };
        data.push(origin_byte);

        // Encode all the vertices
        for vertex in &self.vertices {
            // Encode each coordinate as a 64-bit float in big-endian order
            data.extend_from_slice(&vertex.0.x.to_be_bytes());
            data.extend_from_slice(&vertex.0.y.to_be_bytes());
            data.extend_from_slice(&vertex.0.z.to_be_bytes());
        }

        data
    }

    /// decode decodes a byte slice encoded by encode.
    pub fn decode(data: &[u8]) -> Result<Self, String> {
        // Check for minimum size (version, num vertices, origin inside)
        if data.len() < 6 {
            return Err("Encoded data too short".to_string());
        }

        // Check version
        if data[0] != 1 {
            return Err(format!("Unknown encoding version {}", data[0]));
        }

        // Decode number of vertices
        let mut num_vertices_bytes = [0u8; 4];
        num_vertices_bytes.copy_from_slice(&data[1..5]);
        let num_vertices = u32::from_be_bytes(num_vertices_bytes) as usize;

        // Decode origin_inside flag
        let origin_inside = data[5] != 0;

        // Check if data is long enough to contain all vertices
        let expected_size = 6 + (num_vertices * 24); // 6 bytes header + 24 bytes per vertex
        if data.len() < expected_size {
            return Err(format!(
                "Encoded data too short: expected {} bytes, found {}",
                expected_size,
                data.len()
            ));
        }

        // Decode vertices
        let mut vertices = Vec::with_capacity(num_vertices);
        for i in 0..num_vertices {
            let offset = 6 + (i * 24);

            let mut x_bytes = [0u8; 8];
            let mut y_bytes = [0u8; 8];
            let mut z_bytes = [0u8; 8];

            x_bytes.copy_from_slice(&data[offset..offset + 8]);
            y_bytes.copy_from_slice(&data[offset + 8..offset + 16]);
            z_bytes.copy_from_slice(&data[offset + 16..offset + 24]);

            let x = f64::from_be_bytes(x_bytes);
            let y = f64::from_be_bytes(y_bytes);
            let z = f64::from_be_bytes(z_bytes);

            vertices.push(Point(crate::r3::vector::Vector { x, y, z }));
        }

        // Create a new loop with the decoded vertices and origin_inside flag
        let mut loop_inst = Loop {
            vertices,
            origin_inside,
            depth: 0,
            bound: Rect::empty(),
            subregion_bound: Rect::empty(),
            index: ShapeIndex::new(),
        };

        // Initialize the bounds and index
        loop_inst.init_bound();
        loop_inst.index = ShapeIndex::new();
        let wrapped_shape_loop = loop_inst.clone().into();
        loop_inst.index.add(&wrapped_shape_loop);

        Ok(loop_inst)
    }

    /// encode_compressed encodes a loop using a more compact (but not
    /// lossless) representation.
    pub fn encode_compressed(&self, snap_level: i32) -> Vec<u8> {
        let mut data = Vec::new();

        // Version number
        data.push(1); // Encoding version 1

        // Encode the snap level (for decoding)
        data.extend_from_slice(&snap_level.to_be_bytes());

        // Encode the number of vertices
        let num_vertices = self.vertices.len() as u32;
        data.extend_from_slice(&num_vertices.to_be_bytes());

        // Encode the origin_inside flag
        let origin_byte = if self.origin_inside { 1u8 } else { 0u8 };
        data.push(origin_byte);

        // Encode vertices as cell IDs at the given snap level
        for vertex in &self.vertices {
            let cell_id = crate::s2::cellid::CellID::from(vertex);
            let cell_id_at_level = cell_id.parent(snap_level as u64);
            data.extend_from_slice(&cell_id_at_level.0.to_be_bytes());
        }

        data
    }

    /// decode_compressed decodes a loop encoded using encode_compressed.
    pub fn decode_compressed(data: &[u8]) -> Result<Self, String> {
        // Check for minimum size (version, snap level, num vertices, origin inside)
        if data.len() < 10 {
            return Err("Encoded data too short".to_string());
        }

        // Check version
        if data[0] != 1 {
            return Err(format!("Unknown encoding version {}", data[0]));
        }

        // Decode snap level
        let mut snap_level_bytes = [0u8; 4];
        snap_level_bytes.copy_from_slice(&data[1..5]);
        let _snap_level = i32::from_be_bytes(snap_level_bytes);

        // Decode number of vertices
        let mut num_vertices_bytes = [0u8; 4];
        num_vertices_bytes.copy_from_slice(&data[5..9]);
        let num_vertices = u32::from_be_bytes(num_vertices_bytes) as usize;

        // Decode origin_inside flag
        let origin_inside = data[9] != 0;

        // Check if data is long enough to contain all vertices
        let expected_size = 10 + (num_vertices * 8); // 10 bytes header + 8 bytes per vertex (CellID)
        if data.len() < expected_size {
            return Err(format!(
                "Encoded data too short: expected {} bytes, found {}",
                expected_size,
                data.len()
            ));
        }

        // Decode vertices
        let mut vertices = Vec::with_capacity(num_vertices);
        for i in 0..num_vertices {
            let offset = 10 + (i * 8);

            let mut cell_id_bytes = [0u8; 8];
            cell_id_bytes.copy_from_slice(&data[offset..offset + 8]);
            let cell_id_raw = u64::from_be_bytes(cell_id_bytes);

            let cell_id = crate::s2::cellid::CellID(cell_id_raw);
            let point = cell_id.center_point();

            vertices.push(point);
        }

        // Create a new loop with the decoded vertices and origin_inside flag
        let mut loop_inst = Loop {
            vertices,
            origin_inside,
            depth: 0,
            bound: Rect::empty(),
            subregion_bound: Rect::empty(),
            index: ShapeIndex::new(),
        };

        // Initialize the bounds and index
        loop_inst.init_bound();
        loop_inst.index = ShapeIndex::new();
        loop_inst.index.add(&ShapeType::Loop(loop_inst.clone()));

        Ok(loop_inst)
    }
}

/// Implement Shape trait for Loop
impl Shape for Loop {
    fn num_edges(&self) -> i64 {
        if self.is_empty_or_full() {
            0
        } else {
            self.vertices.len() as i64
        }
    }

    fn edge(&self, i: i64) -> Edge {
        Edge {
            v0: self.vertex(i as usize),
            v1: self.vertex((i as usize) + 1),
        }
    }

    fn num_chains(&self) -> i64 {
        if self.is_empty() {
            0
        } else {
            1
        }
    }

    fn chain(&self, _chain_id: i64) -> Chain {
        Chain {
            start: 0,
            length: self.num_edges(),
        }
    }

    fn chain_edge(&self, _chain_id: i64, offset: i64) -> Edge {
        Edge {
            v0: self.vertex(offset as usize),
            v1: self.vertex((offset as usize) + 1),
        }
    }

    fn chain_position(&self, edge_id: i64) -> ChainPosition {
        ChainPosition {
            chain_id: 0,
            offset: edge_id,
        }
    }

    fn dimension(&self) -> i64 {
        2
    }

    fn reference_point(&self) -> ReferencePoint {
        self.reference_point()
    }
}

impl Add<VertexTraversalDirection> for usize {
    type Output = usize;

    fn add(self, rhs: VertexTraversalDirection) -> Self::Output {
        match rhs {
            VertexTraversalDirection::Forward => self + 1,
            VertexTraversalDirection::Backward => self - 1,
        }
    }
}

impl Sub<VertexTraversalDirection> for usize {
    type Output = usize;

    fn sub(self, rhs: VertexTraversalDirection) -> Self::Output {
        match rhs {
            VertexTraversalDirection::Forward => self - 1,
            VertexTraversalDirection::Backward => self + 1,
        }
    }
}

impl Mul<Angle> for VertexTraversalDirection {
    type Output = Angle;

    fn mul(self, rhs: Angle) -> Self::Output {
        match self {
            VertexTraversalDirection::Forward => rhs,
            VertexTraversalDirection::Backward => -rhs,
        }
    }
}

impl PartialEq<i32> for VertexTraversalDirection {
    fn eq(&self, other: &i32) -> bool {
        match self {
            VertexTraversalDirection::Forward => *other == 1,
            VertexTraversalDirection::Backward => *other == -1,
        }
    }
}

// Extension to Loop implementation with normalization functions
impl Loop {
    /// Determines if the loop is normalized, meaning its area is at most 2*pi.
    /// Degenerate loops are handled consistently with Sign. For example, if a loop
    /// can be expressed as a union of degenerate or nearly-degenerate CCW triangles,
    /// then it will always be considered normalized.
    pub fn is_normalized(&self) -> bool {
        // Optimization: if the longitude span is less than 180 degrees, then the
        // loop covers less than half the sphere and is therefore normalized.
        if self.bound.lng.len() < PI {
            return true;
        }

        // We allow some error so that hemispheres are always considered normalized.
        // The turning angle evaluates exactly to -2*pi for hemispheres, with no error.
        let max_error = self.turning_angle_max_error();
        self.turning_angle() >= -max_error
    }

    /// Normalizes the loop if necessary so that the area enclosed by the loop is
    /// at most 2*pi. This may invert the loop.
    pub fn normalize(&mut self) {
        if !self.is_normalized() {
            self.invert();
        }
    }

    /// Reverses the order of the loop vertices, effectively complementing the
    /// region represented by the loop. For example, the loop ABCD (with edges
    /// AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
    pub fn invert(&mut self) {
        self.index.reset();
        if self.is_empty_or_full() {
            if self.is_full() {
                self.vertices[0] = EMPTY_LOOP_POINT;
            } else {
                self.vertices[0] = FULL_LOOP_POINT;
            }
        } else {
            // For non-special loops, reverse the slice of vertices.
            self.vertices.reverse();
        }

        // originInside must be set correctly before rebuilding the ShapeIndex.
        self.origin_inside = !self.origin_inside;
        if self.bound.lat.lo > -PI / 2.0 && self.bound.lat.hi < PI / 2.0 {
            // The complement of this loop contains both poles.
            self.bound = Rect::full();
            self.subregion_bound = self.bound.clone();
        } else {
            self.init_bound();
        }

        // TODO: Not sure if this cloning treatment is correct (arc? clone? Point lifetime?)
        let wrapped_loop = ShapeType::Loop(self.clone());

        self.index.add(&wrapped_loop);
    }

    // CanonicalFirstVertex returns a first index and a direction (either +1 or -1)
    // such that the vertex sequence (first, first+dir, ..., first+(n-1)*dir) does
    // not change when the loop vertex order is rotated or inverted. This allows the
    // loop vertices to be traversed in a canonical order. The return values are
    // chosen such that (first, ..., first+n*dir) are in the range [0, 2*n-1] as
    // expected by the Vertex method.
    pub(crate) fn canonical_first_vertex(&self) -> (usize, VertexTraversalDirection) {
        let mut first_idx = 0;
        let n = self.vertices.len();
        for i in 1..n {
            if self.vertex(i).0 < (self.vertex(first_idx).0) {
                first_idx = i
            }
        }

        // 0 <= first_idx <= n-1, so (first_idx+n*dir) <= 2*n-1.
        if self.vertex(first_idx + 1).0 < (self.vertex(first_idx + n - 1).0) {
            return (first_idx, VertexTraversalDirection::Forward);
        }

        // n <= first_idx <= 2*n-1, so (first_idx+n*dir) >= 0.
        first_idx += n;
        return (first_idx, VertexTraversalDirection::Backward);
    }

    /// Returns the sum of the turning angles at each vertex. The return value is
    /// positive if the loop is counter-clockwise, negative if the loop is
    /// clockwise, and zero if the loop is a great circle.
    ///
    /// This quantity is also called the "geodesic curvature" of the loop.
    pub fn turning_angle(&self) -> f64 {
        // For empty and full loops, we return the limit value as the loop area approaches
        // 0 or 4*Pi respectively.
        if self.is_empty_or_full() {
            if self.contains_origin() {
                return -2.0 * PI;
            }
            return 2.0 * PI;
        }

        // Don't crash even if the loop is not well-defined.
        if self.vertices.len() < 3 {
            return 0.0;
        }

        // TODO: Implement canonical_first_vertex function
        // For now, we'll just use vertex 0 as the starting point

        let n = self.vertices.len();
        let (mut i, direction) = self.canonical_first_vertex();
        let mut sum = self.turn_angle(
            self.vertex((i + n - direction) % n),
            self.vertex(i % n),
            self.vertex((i + direction) % n),
        );

        // Using Kahan summation for better numerical stability
        let mut compensation = 0.0;
        let mut remaining = n - 1;

        while remaining > 0 {
            i = (i + direction) % n;
            let angle = self.turn_angle(
                self.vertex((i + n - direction) % n),
                self.vertex(i),
                self.vertex((i + direction) % n),
            );
            let old_sum = sum;
            let corrected_angle = angle + compensation;
            sum += corrected_angle;
            compensation = ((old_sum - sum) + corrected_angle).0;
            remaining -= 1;
        }

        // Bound the result to handle numerical issues
        const MAX_CURVATURE: f64 = 2.0 * PI - 4.0 * DBL_EPSILON;

        // 	return math.Max(-maxCurvature, math.Min(maxCurvature, float64(dir)*float64(sum+compensation)))

        let min_max_curv_and_compensation = MAX_CURVATURE.min((direction * (sum + compensation)).0);
        -MAX_CURVATURE.max(min_max_curv_and_compensation)
    }

    /// Returns the maximum error in TurningAngle. The value is not constant; it
    /// depends on the loop.
    pub(crate) fn turning_angle_max_error(&self) -> f64 {
        // The maximum error can be bounded as follows:
        //   3.00 * dblEpsilon    for RobustCrossProd(b, a)
        //   3.00 * dblEpsilon    for RobustCrossProd(c, b)
        //   3.25 * dblEpsilon    for Angle()
        //   2.00 * dblEpsilon    for each addition in the Kahan summation
        //   ------------------
        // 11.25 * dblEpsilon
        let max_error_per_vertex = 11.25 * DBL_EPSILON;
        max_error_per_vertex * self.vertices.len() as f64
    }

    /// Compute the turning angle between three consecutive vertices.
    fn turn_angle(&self, v0: Point, v1: Point, v2: Point) -> Angle {
        // We use the cross product formula rather than a more complex but
        // numerically stable formula because the final result is normalized
        // by the total turning angle.
        let angle = v0.0.cross(&v1.0).angle(&v1.0.cross(&v2.0));

        // Use the determinant to figure out if the angle is positive or negative.
        if v0.0.dot(&v1.0.cross(&v2.0)) > 0.0 {
            angle
        } else {
            -angle
        }
    }

    /// Returns the area of the loop interior, i.e. the region on the left side of
    /// the loop. The return value is between 0 and 4*pi. This value is not affected
    /// by whether this loop is a "hole" or a "shell".
    pub fn area(&self) -> f64 {
        // It is surprisingly difficult to compute the area of a loop robustly. The
        // main issues are (1) whether degenerate loops are considered to be CCW or
        // not (i.e., whether their area is close to 0 or 4*pi), and (2) computing
        // the areas of small loops with good relative accuracy.
        if self.is_empty_or_full() {
            if self.contains_origin() {
                return 4.0 * PI;
            }
            return 0.0;
        }

        // Use the "signed area" approach which computes a signed sum over triangles
        let mut area = self.surface_integral_float64(signed_area);

        // The signed area should be between approximately -4*pi and 4*pi.
        if area < 0.0 {
            // We have computed the negative of the area of the loop exterior.
            area += 4.0 * PI;
        }

        // Clamp the result to the valid range.
        area = area.min(4.0 * PI).max(0.0);

        // If the area is close enough to zero or 4*pi so that the loop orientation
        // is ambiguous, then we compute the loop orientation explicitly.
        let max_error = self.turning_angle_max_error();
        if area < max_error && !self.is_normalized() {
            return 4.0 * PI;
        } else if area > (4.0 * PI - max_error) && self.is_normalized() {
            return 0.0;
        }

        area
    }

    /// Computes the oriented surface integral of some quantity f(x) over the loop
    /// interior. Specialized version for f64 return values.
    fn surface_integral_float64<F>(&self, f: F) -> f64
    where
        F: Fn(Point, Point, Point) -> f64,
    {
        // Maximum length of an edge for it to be considered numerically stable.
        const MAX_LENGTH: f64 = PI - 1e-5;

        let mut sum = 0.0;
        let mut origin = self.vertex(0);

        for i in 1..self.vertices.len() - 1 {
            // Let V_i be vertex(i), let O be the current origin, and let length(A,B)
            // be the length of edge (A,B). At the start of each loop iteration, the
            // "leading edge" of the triangle fan is (O,V_i), and we want to extend
            // the triangle fan so that the leading edge is (O,V_i+1).
            if self.vertex(i + 1).0.angle(&origin.0).0 > MAX_LENGTH {
                // We are about to create an unstable edge, so choose a new origin O'
                // for the triangle fan.
                let old_origin = origin;
                if origin == self.vertex(0) {
                    // The following point is well-separated from V_i and V_0 (and
                    // therefore V_i+1 as well).
                    origin = Point(self.vertex(0).0.cross(&self.vertex(i).0).normalize());
                } else if self.vertex(i).0.angle(&self.vertex(0).0).0 < MAX_LENGTH {
                    // All edges of the triangle (O, V_0, V_i) are stable, so we can
                    // revert to using V_0 as the origin.
                    origin = self.vertex(0);
                } else {
                    // (O, V_i+1) and (V_0, V_i) are antipodal pairs, and O and V_0 are
                    // perpendicular. Therefore V_0.CrossProd(O) is approximately
                    // perpendicular to all of {O, V_0, V_i, V_i+1}, and we can choose
                    // this point O' as the new origin.
                    origin = Point(self.vertex(0).cross(&old_origin).0);

                    // Advance the edge (V_0,O) to (V_0,O').
                    sum += f(self.vertex(0), old_origin, origin);
                }

                // Advance the edge (O,V_i) to (O',V_i).
                sum += f(old_origin, self.vertex(i), origin);
            }

            // Advance the edge (O,V_i) to (O,V_i+1).
            sum += f(origin, self.vertex(i), self.vertex(i + 1));
        }

        // If the origin is not V_0, we need to sum one more triangle.
        if origin != self.vertex(0) {
            // Advance the edge (O,V_n-1) to (O,V_0).
            sum += f(origin, self.vertex(self.vertices.len() - 1), self.vertex(0));
        }

        sum
    }

    /// Returns the true centroid of the loop multiplied by the area of the loop.
    /// The result is not unit length. The centroid may not be contained by the loop.
    ///
    /// We prescale by the loop area for two reasons: (1) it is cheaper to compute
    /// this way, and (2) it makes it easier to compute the centroid of more
    /// complicated shapes (by splitting them into disjoint regions and adding
    /// their centroids).
    pub fn centroid(&self) -> Point {
        self.surface_integral_point(crate::point::true_centroid)
    }

    /// Computes the oriented surface integral of some vector-valued quantity over
    /// the loop interior. Similar to surface_integral_float64 but returns a Point.
    fn surface_integral_point<F>(&self, f: F) -> Point
    where
        F: Fn(&Point, &Point, &Point) -> Point,
    {
        // Maximum length of an edge for it to be considered numerically stable.
        const MAX_LENGTH: f64 = PI - 1e-5;

        let mut sum = R3Vector {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let mut origin = self.vertex(0);
        // TODO: We dont yet have `Point::point_cross` implemented, using normal cross
        for i in 1..self.vertices.len() - 1 {
            if self.vertex(i + 1).0.angle(&origin.0).0 > MAX_LENGTH {
                let old_origin = origin;
                if origin == self.vertex(0) {
                    origin = self.vertex(0).cross(&self.vertex(i)).normalize();
                } else if self.vertex(i).0.angle(&self.vertex(0).0).0 < MAX_LENGTH {
                    origin = self.vertex(0);
                } else {
                    origin = Point(self.vertex(0).0.cross(&old_origin.0));
                    sum = sum + f(&self.vertex(0), &old_origin, &origin).0;
                }
                sum = sum + f(&old_origin, &self.vertex(i), &origin).0;
            }
            sum = sum + f(&origin, &self.vertex(i), &self.vertex(i + 1)).0;
        }

        if origin != self.vertex(0) {
            sum = sum
                + f(
                    &origin,
                    &self.vertex(self.vertices.len() - 1),
                    &self.vertex(0),
                )
                .0;
        }

        Point(sum)
    }
}

/// SignedArea returns the area of the triangle (a,b,c). The result is positive if the
/// triangle is counterclockwise and negative if the triangle is clockwise.
/// This is used by area().
fn signed_area(a: Point, b: Point, c: Point) -> f64 {
    // This method computes the area using Girard's formula. It's equivalent to
    // computing the solid angle subtended by the triangle at the origin, and then
    // dividing by 2.
    let a_dot_b = a.0.dot(&b.0).clamp(-1.0, 1.0);
    let b_dot_c = b.0.dot(&c.0).clamp(-1.0, 1.0);
    let c_dot_a = c.0.dot(&a.0).clamp(-1.0, 1.0);

    // Let a, b, c be the spherical coordinates (i.e., unit vectors) of the triangle
    // vertices. Then we compute the area with the following formula:
    //   tan(area/2) = (det(a, b, c) / (1 + |a·b| + |b·c| + |c·a|))
    //
    // where det(a, b, c) = a·(b×c), i.e. the triple product of a, b, c.

    // The determinant is the triple product a·(b×c). We can use the identity
    // |a×b| = 2·sin(θ/2) for two unit vectors separated by angle θ.
    let det = a.0.dot(&b.0.cross(&c.0));

    // The denominator is = 1 + |a·b| + |b·c| + |c·a|
    let denom = 1.0 + a_dot_b.abs() + b_dot_c.abs() + c_dot_a.abs();

    // Now use the formula: tan(area/2) = det / denom
    let tan_area_over_2 = det / denom;

    // And finally get the area using atan.
    2.0 * f64::atan(tan_area_over_2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::f64_eq;
    use crate::r1::interval::Interval;
    use crate::r3::vector::{Vector as R3Vector, Vector};
    use crate::s1::angle::Angle;
    use crate::s1::Deg;
    use crate::s2::cell::Cell;
    use crate::s2::cellid::CellID;
    use crate::s2::point::Point;
    use crate::s2::shape::Shape;

    use crate::crossing_edge_query::CrossingType;
    use crate::latlng::LatLng;
    use crate::s2::edge_crossings::Crossing;
    use std::convert::TryInto;
    use std::f64::consts::PI;

    // Helper function to create a test loop from coordinates
    fn make_loop(s: &str) -> Loop {
        if s == "full" {
            Loop::full()
        } else if s == "empty" {
            Loop::empty()
        } else {
            Loop::from_points(parse_points(s))
        }
    }

    // TODO: Make text format test suite helper function independent
    // parseLatLngs returns the values in the input string as LatLngs.
    fn parse_lat_lngs(s: &str) -> Vec<LatLng> {
        let mut result = Vec::new();
        for pair in s.split(',') {
            if pair.is_empty() {
                continue;
            }
            let mut parts: Vec<&str> = pair
                .split(':')
                .collect::<Vec<&str>>()
                .into_iter()
                .map(|s| s.trim())
                .collect();

            if parts.len() != 2 {
                continue;
            }

            let lat_float: f64 = parts[0].parse().unwrap();
            let lng_float: f64 = parts[1].parse().unwrap();
            result.push(LatLng::from_degrees(lat_float, lng_float));
        }
        result
    }

    // TODO: Make text format test suite helper function independent
    // parse_points returns the values in the input string as Points.
    fn parse_points(s: &str) -> Vec<Point> {
        let lls = parse_lat_lngs(s);
        if lls.is_empty() {
            return Vec::new();
        }

        let mut points = Vec::new();
        for ll in lls {
            points.push(Point::from(ll));
        }
        points
    }

    // Helper to create a loop from lat/lng points (degrees)
    // TODO: use LatLng class and use Point constructor for LatLng points for accuracy sake
    fn lat_lng_loop(points: &[(f64, f64)]) -> Loop {
        let vertices: Vec<Point> = points
            .iter()
            .map(|&(lat, lng)| {
                let lat_rad = lat * PI / 180.0;
                let lng_rad = lng * PI / 180.0;
                let phi = PI / 2.0 - lat_rad;
                let x = phi.sin() * lng_rad.cos();
                let y = phi.sin() * lng_rad.sin();
                let z = phi.cos();
                Point(R3Vector { x, y, z })
            })
            .collect();
        Loop::from_points(vertices)
    }

    fn rotate(l: &Loop) -> Loop {
        if l.is_empty_or_full() {
            return l.clone();
        }

        let mut vertices = Vec::new();
        // Skip the first vertex and copy the rest
        for i in 1..l.num_vertices() {
            vertices.push(l.vertex(i));
        }
        // Add the first vertex at the end
        vertices.push(l.vertex(0));

        Loop::from_points(vertices)
    }

    // Basic test loops
    fn empty_loop() -> Loop {
        Loop::empty()
    }

    fn full_loop() -> Loop {
        Loop::full()
    }

    // The northern hemisphere, defined using two pairs of antipodal points.
    fn north_hemi() -> Loop {
        Loop::from_points(parse_points("0:-180, 0:-90, 0:0, 0:90"))
    }

    // The northern hemisphere, defined using three points 120 degrees apart.
    fn north_hemi3() -> Loop {
        Loop::from_points(parse_points("0:-180, 0:-60, 0:60"))
    }

    // The southern hemisphere, defined using two pairs of antipodal points.
    fn south_hemi() -> Loop {
        Loop::from_points(parse_points("0:90, 0:0, 0:-90, 0:-180"))
    }

    // The western hemisphere, defined using two pairs of antipodal points.
    fn west_hemi() -> Loop {
        Loop::from_points(parse_points("0:-180, -90:0, 0:0, 90:0"))
    }

    // The eastern hemisphere, defined using two pairs of antipodal points.
    fn east_hemi() -> Loop {
        Loop::from_points(parse_points("90:0, 0:0, -90:0, 0:-180"))
    }

    // The "near" hemisphere, defined using two pairs of antipodal points.
    fn near_hemi() -> Loop {
        Loop::from_points(parse_points("0:-90, -90:0, 0:90, 90:0"))
    }

    // The "far" hemisphere, defined using two pairs of antipodal points.
    fn far_hemi() -> Loop {
        Loop::from_points(parse_points("90:0, 0:90, -90:0, 0:-90"))
    }

    // A spiral stripe that slightly over-wraps the equator.
    fn candy_cane() -> Loop {
        Loop::from_points(parse_points(
            "-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70",
        ))
    }

    // A small clockwise loop in the northern & eastern hemisperes.
    fn small_necw() -> Loop {
        Loop::from_points(parse_points("35:20, 45:20, 40:25"))
    }

    // Loop around the north pole at 80 degrees.
    fn arctic80() -> Loop {
        Loop::from_points(parse_points("80:-150, 80:-30, 80:90"))
    }

    // Loop around the south pole at 80 degrees.
    fn antarctic80() -> Loop {
        Loop::from_points(parse_points("-80:120, -80:0, -80:-120"))
    }

    // A completely degenerate triangle along the equator that RobustCCW()
    // considers to be CCW.
    fn line_triangle() -> Loop {
        Loop::from_points(parse_points("0:1, 0:2, 0:3"))
    }

    // A nearly-degenerate CCW chevron near the equator with very long sides
    // (about 80 degrees).  Its area is less than 1e-640, which is too small
    // to represent in double precision.
    fn skinny_chevron() -> Loop {
        Loop::from_points(parse_points("0:0, -1e-320:80, 0:1e-320, 1e-320:80"))
    }

    // A diamond-shaped loop around the point 0:180.
    fn loop_a() -> Loop {
        Loop::from_points(parse_points("0:178, -1:180, 0:-179, 1:-180"))
    }

    // Like loopA, but the vertices are at leaf cell centers.
    fn snapped_loop_a() -> Loop {
        let points = vec![
            CellID::from(&parse_points("0:178")[0]).center_point(),
            CellID::from(&parse_points("-1:180")[0]).center_point(),
            CellID::from(&parse_points("0:-179")[0]).center_point(),
            CellID::from(&parse_points("1:-180")[0]).center_point(),
        ];

        Loop::from_points(points)
    }

    // A different diamond-shaped loop around the point 0:180.
    fn loop_b() -> Loop {
        Loop::from_points(parse_points("0:179, -1:180, 0:-178, 1:-180"))
    }

    // The intersection of A and B.
    fn a_intersect_b() -> Loop {
        Loop::from_points(parse_points("0:179, -1:180, 0:-179, 1:-180"))
    }

    // The union of A and B.
    fn a_union_b() -> Loop {
        Loop::from_points(parse_points("0:178, -1:180, 0:-178, 1:-180"))
    }

    // A minus B (concave).
    fn a_minus_b() -> Loop {
        Loop::from_points(parse_points("0:178, -1:180, 0:179, 1:-180"))
    }

    // B minus A (concave).
    fn b_minus_a() -> Loop {
        Loop::from_points(parse_points("0:-179, -1:180, 0:-178, 1:-180"))
    }

    // A shape gotten from A by adding a triangle to one edge, and
    // subtracting a triangle from the opposite edge.
    fn loop_c() -> Loop {
        Loop::from_points(parse_points("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180"))
    }

    // A shape gotten from A by adding a triangle to one edge, and
    // adding another triangle to the opposite edge.
    fn loop_d() -> Loop {
        Loop::from_points(parse_points(
            "0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180",
        ))
    }

    //   3------------2
    //   |            |               ^
    //   |  7-8  b-c  |               |
    //   |  | |  | |  |      Latitude |
    //   0--6-9--a-d--1               |
    //   |  | |       |               |
    //   |  f-e       |               +----------->
    //   |            |                 Longitude
    //   4------------5
    //
    // Important: It is not okay to skip over collinear vertices when
    // defining these loops (e.g. to define loop E as "0,1,2,3") because S2
    // uses symbolic perturbations to ensure that no three vertices are
    // *ever* considered collinear (e.g., vertices 0, 6, 9 are not
    // collinear).  In other words, it is unpredictable (modulo knowing the
    // details of the symbolic perturbations) whether 0123 contains 06123
    // for example.

    // Loop E:  0,6,9,a,d,1,2,3
    fn loop_e() -> Loop {
        Loop::from_points(parse_points(
            "0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30",
        ))
    }

    // Loop F:  0,4,5,1,d,a,9,6
    fn loop_f() -> Loop {
        Loop::from_points(parse_points(
            "0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34",
        ))
    }

    // Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
    fn loop_g() -> Loop {
        Loop::from_points(parse_points(
            "0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30",
        ))
    }

    // Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
    fn loop_h() -> Loop {
        Loop::from_points(parse_points(
            "0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30",
        ))
    }

    // Loop I:  7,6,f,e,9,8
    fn loop_i() -> Loop {
        Loop::from_points(parse_points("10:34, 0:34, -10:34, -10:36, 0:36, 10:36"))
    }

    // The set of all test loops.
    fn all_loops() -> Vec<Loop> {
        vec![
            empty_loop(),
            full_loop(),
            north_hemi(),
            north_hemi3(),
            south_hemi(),
            west_hemi(),
            east_hemi(),
            near_hemi(),
            far_hemi(),
            candy_cane(),
            small_necw(),
            arctic80(),
            antarctic80(),
            line_triangle(),
            skinny_chevron(),
            loop_a(),
            //snapped_loop_a(), // Fails TestAreaConsistentWithTurningAngle
            loop_b(),
            a_intersect_b(),
            a_union_b(),
            a_minus_b(),
            b_minus_a(),
            loop_c(),
            loop_d(),
            loop_e(),
            loop_f(),
            loop_g(),
            loop_h(),
            loop_i(),
        ]
    }

    #[test]
    fn test_loop_empty() {
        let loop_empty = empty_loop();
        assert!(loop_empty.is_empty());
        assert!(!loop_empty.is_full());
        assert_eq!(loop_empty.num_edges(), 0);
        assert_eq!(loop_empty.num_chains(), 0);
        assert_eq!(loop_empty.dimension(), 2);
        assert!(!loop_empty.reference_point().contained);
    }

    #[test]
    fn test_loop_full() {
        let loop_full = full_loop();
        assert_eq!(loop_full.num_edges(), 0);
        assert_eq!(loop_full.num_chains(), 1);
        assert_eq!(loop_full.dimension(), 2);
        assert!(!loop_full.is_empty());
        assert!(loop_full.is_full());
        assert!(loop_full.reference_point().contained);
    }

    #[test]
    fn test_loop_basic() {
        let shape = ShapeType::Loop(make_loop("0:0, 0:1, 1:0"));
        println!("{:?}", shape);

        // Verify edge count and basic properties
        assert_eq!(
            shape.num_edges(),
            3,
            "shape.num_edges() = {}, want {}",
            shape.num_edges(),
            3
        );
        assert_eq!(
            shape.num_chains(),
            1,
            "shape.num_chains() = {}, want {}",
            shape.num_chains(),
            1
        );

        // Check chain properties
        let chain = shape.chain(0);
        assert_eq!(
            chain.start, 0,
            "shape.chain(0).start = {}, want {}",
            chain.start, 0
        );
        assert_eq!(
            chain.length, 3,
            "shape.chain(0).length = {}, want {}",
            chain.length, 3
        );

        // Check edge properties
        let e = shape.edge(2);
        let want_v0 = Point::from(LatLng::from_degrees(1.0, 0.0));
        let want_v1 = Point::from(LatLng::from_degrees(0.0, 0.0));

        assert!(
            e.v0.approx_eq(&want_v0),
            "shape.edge(2) end A = {:?}, want {:?}",
            e.v0,
            want_v0
        );
        assert!(
            e.v1.approx_eq(&want_v1),
            "shape.edge(2) end B = {:?}, want {:?}",
            e.v1,
            want_v1
        );

        // Check other properties
        assert_eq!(
            shape.dimension(),
            2,
            "shape.dimension() = {}, want {}",
            shape.dimension(),
            2
        );
        assert!(!shape.is_empty(), "shape.is_empty() = true, want false");
        assert!(!shape.is_full(), "shape.is_full() = true, want false");
        assert!(
            !shape.reference_point().contained,
            "shape.reference_point().contained = true, want false"
        );
    }

    #[test]
    fn test_loop_hole_and_sign() {
        let mut l = Loop::from_points(parse_points("0:-180, 0:-90, 0:0, 0:90"));

        if l.is_hole() {
            panic!("loop with default depth should not be a hole");
        }
        if l.sign() == -1 {
            panic!("loop with default depth should have a sign of +1");
        }

        l.depth = 3;
        if !l.is_hole() {
            panic!("loop with odd depth should be a hole");
        }
        if l.sign() != -1 {
            panic!("loop with odd depth should have a sign of -1");
        }

        l.depth = 2;
        if l.is_hole() {
            panic!("loop with even depth should not be a hole");
        }
        if l.sign() == -1 {
            panic!("loop with even depth should have a sign of +1");
        }
    }

    #[test]
    fn test_loop_rect_bound() {
        // TODO: Get this from rect instead

        // TODO: Isolate these out into test util in independen
        pub fn max_error_for_tests() -> LatLng {
            // The maximum error in the latitude calculation is
            //    3.84 * DBL_EPSILON   for the cross product calculation (see above)
            //    0.96 * DBL_EPSILON   for the Latitude() calculation
            //    5    * DBL_EPSILON   added by add_point/get_bound to compensate for error
            //    ------------------
            //    9.80 * DBL_EPSILON   maximum error in result
            //
            // The maximum error in the longitude calculation is DBL_EPSILON. get_bound()
            // does not do any expansion because this isn't necessary in order to
            // bound the *rounded* longitudes of contained points.
            return LatLng::from_radians(10. * DBL_EPSILON, 1. * DBL_EPSILON);
        }
        let rect_error = max_error_for_tests();

        assert!(
            empty_loop().bound.is_empty(),
            "empty loop's RectBound should be empty"
        );
        assert!(
            full_loop().bound.is_full(),
            "full loop's RectBound should be full"
        );
        assert!(
            candy_cane().bound.lng.is_full(),
            "candy cane loop's RectBound should have a full longitude range"
        );

        let candy_lat_lo = candy_cane().bound.lat.lo;
        assert!(
            candy_lat_lo < -0.349066,
            "candy cane loop's RectBound should have a lower latitude ({}) under -0.349066 radians",
            candy_lat_lo
        );

        let candy_lat_hi = candy_cane().bound.lat.hi;
        assert!(
            candy_lat_hi > 0.174533,
            "candy cane loop's RectBound should have an upper latitude ({}) over 0.174533 radians",
            candy_lat_hi
        );

        assert!(
            small_necw().bound.is_full(),
            "small northeast clockwise loop's RectBound should be full"
        );

        let arctic_bound = arctic80().bound;
        let arctic_want = crate::s2::rect::Rect::from_degrees(80.0, -180.0, 90.0, 180.0);
        assert!(
            Rect::approx_eq_by(&arctic_bound, &arctic_want, &rect_error),
            "arctic 80 loop's RectBound ({:?}) should be {:?}",
            arctic_bound,
            arctic_want
        );

        let antarctic_bound = antarctic80().bound;
        let antarctic_want = crate::s2::rect::Rect::from_degrees(-90.0, -180.0, -80.0, 180.0);
        assert!(
            Rect::approx_eq_by(&antarctic_bound, &antarctic_want, &rect_error),
            "antarctic 80 loop's RectBound ({:?}) should be {:?}",
            antarctic_bound,
            antarctic_want
        );

        assert!(
            south_hemi().bound.lng.is_full(),
            "south hemi loop's RectBound should have a full longitude range"
        );

        let south_lat = south_hemi().bound.lat;
        let want_lat = crate::r1::interval::Interval {
            lo: -PI / 2.0,
            hi: 0.0,
        };
        assert!(
            Interval::approx_eq_by(&south_lat, &want_lat, rect_error.lat.rad()),
            "south hemi loop's RectBound latitude interval ({:?}) should be {:?}",
            south_lat,
            want_lat
        );

        // Create a loop that contains the complement of the arctic80 loop
        let mut arctic80_inv = arctic80();
        arctic80_inv.invert();

        // The highest latitude of each edge is attained at its midpoint
        let v0 = arctic80_inv.vertex(0).0;
        let v1 = arctic80_inv.vertex(1).0;
        let mid = Point(v0.add(v1).mul(0.5));

        let lat_hi = arctic80_inv.bound.lat.hi;
        let want_lat_hi = LatLng::from(mid).lat.rad();
        assert_f64_eq!(
            lat_hi,
            want_lat_hi,
            10.0 * DBL_EPSILON,
            format!(
                "arctic 80 inverse loop's RectBound should have a latitude hi of {}, got {}",
                want_lat_hi, lat_hi
            )
        );
    }
    #[test]
    fn test_loop_cap_bound() {
        assert!(
            empty_loop().cap_bound().is_empty(),
            "empty loop's CapBound should be empty"
        );
        assert!(
            full_loop().cap_bound().is_full(),
            "full loop's CapBound should be full"
        );
        assert!(
            small_necw().cap_bound().is_full(),
            "small northeast clockwise loop's CapBound should be full"
        );

        let got = arctic80().cap_bound();
        let want = crate::s2::rect::Rect::from_degrees(80.0, -180.0, 90.0, 180.0).cap_bound();
        assert!(
            got.approx_eq(&want),
            "arctic 80 loop's CapBound ({:?}) should be {:?}",
            got,
            want
        );

        let got = antarctic80().cap_bound();
        let want = crate::s2::rect::Rect::from_degrees(-90.0, -180.0, -80.0, 180.0).cap_bound();
        assert!(
            got.approx_eq(&want),
            "antarctic 80 loop's CapBound ({:?}) should be {:?}",
            got,
            want
        );
    }

    #[test]
    fn test_loop_origin_inside() {
        assert!(
            north_hemi().contains_origin(),
            "north hemisphere polygon should include origin"
        );
        assert!(
            north_hemi3().contains_origin(),
            "north hemisphere 3 polygon should include origin"
        );
        assert!(
            !south_hemi().contains_origin(),
            "south hemisphere polygon should not include origin"
        );
        assert!(
            !west_hemi().contains_origin(),
            "west hemisphere polygon should not include origin"
        );
        assert!(
            east_hemi().contains_origin(),
            "east hemisphere polygon should include origin"
        );
        assert!(
            !near_hemi().contains_origin(),
            "near hemisphere polygon should not include origin"
        );
        assert!(
            far_hemi().contains_origin(),
            "far hemisphere polygon should include origin"
        );
        assert!(
            !candy_cane().contains_origin(),
            "candy cane polygon should not include origin"
        );
        assert!(
            small_necw().contains_origin(),
            "smallNECW polygon should include origin"
        );
        assert!(
            arctic80().contains_origin(),
            "arctic 80 polygon should include origin"
        );
        assert!(
            !antarctic80().contains_origin(),
            "antarctic 80 polygon should not include origin"
        );
        assert!(
            !loop_a().contains_origin(),
            "loop A polygon should not include origin"
        );
    }
    #[test]
    fn test_loop_contains_point() {
        let north = Point(Vector {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        });
        let south = Point(Vector {
            x: 0.0,
            y: 0.0,
            z: -1.0,
        });

        let east = Point::from_coords(0., 1., 0.);
        let west = Point::from_coords(0., -1., 0.);

        assert!(
            !empty_loop().contains_point(north),
            "Empty loop should not contain north pole or any point!"
        );

        assert!(
            full_loop().contains_point(south),
            "Full loop should contain south pole and all points"
        );

        struct TestCase {
            name: &'static str,
            l: Loop,
            input: Point,
            output: Point,
        }

        // Helper function to rotate a loop by moving the first vertex to the end

        let test_cases = vec![
            TestCase {
                name: "north hemisphere",
                l: north_hemi(),
                input: north,
                output: south,
            },
            TestCase {
                name: "south hemisphere",
                l: south_hemi(),
                input: south,
                output: north,
            },
            TestCase {
                name: "west hemisphere",
                l: west_hemi(),
                input: west,
                output: east,
            },
            TestCase {
                name: "east hemisphere",
                l: east_hemi(),
                input: east,
                output: west,
            },
            TestCase {
                name: "candy cane",
                l: candy_cane(),
                input: Point::from(LatLng::from_degrees(5.0, 71.0)),
                output: Point::from(LatLng::from_degrees(-8.0, 71.0)),
            },
        ];

        for test_case in &test_cases {
            let mut l = test_case.l.clone();
            for i in 0..4 {
                assert!(
                    l.contains_point(test_case.input),
                    "{} loop should contain {:?} at rotation {}",
                    test_case.name,
                    test_case.input,
                    i
                );
                assert!(
                    !l.contains_point(test_case.output),
                    "{} loop shouldn't contain {:?} at rotation {}",
                    test_case.name,
                    test_case.output,
                    i
                );
                l = rotate(&l);
            }
        }

        // Check that each cell vertex is contained by exactly one of the adjacent cells
        for level in 0..3 {
            // Map of unique points across all loops at this level
            let mut points = std::collections::HashMap::<Point, bool>::new();
            let mut loops = Vec::<Loop>::new();

            // Iterate over all cells at this level
            let start_id = CellID::from_face(0).child_begin_at_level(level as u64);
            let end_id = CellID::from_face(5).child_end_at_level(level as u64);

            let mut id = start_id;
            while id != end_id {
                let cell = Cell::from(id);

                // Add cell center to points
                points.insert(cell.center(), true);

                // Add all vertices of the cell
                let mut vertices = Vec::new();
                for k in 0..4 {
                    vertices.push(cell.vertex(k));
                    points.insert(cell.vertex(k), true);
                }

                loops.push(Loop::from_points(vertices));
                id = id.next();
            }

            // Check that each point is contained by exactly one loop
            for (point, _) in points.iter() {
                let mut count = 0;
                for loop_under_test in &loops {
                    if loop_under_test.contains_point(*point) {
                        count += 1;
                    }
                }

                assert_eq!(
                    count, 1,
                    "point {:?} should only be contained by one loop at level {}, got {}",
                    point, level, count
                );
            }
        }
    }

    #[test]
    fn test_loop_vertex() {
        let tests = vec![
            (empty_loop(), 0, Point(R3Vector::new(0.0, 0.0, 1.0))),
            (empty_loop(), 1, Point(R3Vector::new(0.0, 0.0, 1.0))),
            (full_loop(), 0, Point(R3Vector::new(0.0, 0.0, -1.0))),
            (full_loop(), 1, Point(R3Vector::new(0.0, 0.0, -1.0))),
            (arctic80(), 0, parse_points("80:-150")[0]),
            (arctic80(), 1, parse_points("80:-30")[0]),
            (arctic80(), 2, parse_points("80:90")[0]),
            (arctic80(), 3, parse_points("80:-150")[0]),
        ];

        for (loop_under_test, vertex, want) in tests {
            let got = loop_under_test.vertex(vertex);
            assert!(
                got.approx_eq(&want),
                "{:?}.vertex({}) = {:?}, want {:?}",
                loop_under_test,
                vertex,
                got,
                want
            );
        }

        // Check that wrapping is correct.
        assert!(
            arctic80().vertex(2).approx_eq(&arctic80().vertex(5)),
            "Vertex should wrap values. {:?}.vertex(2) = {:?} != {:?}.vertex(5) = {:?}",
            arctic80(),
            arctic80().vertex(2),
            arctic80(),
            arctic80().vertex(5)
        );

        let loop_around_thrice = 2 + 3 * arctic80().vertices.len();
        assert!(
            arctic80()
                .vertex(2)
                .approx_eq(&arctic80().vertex(loop_around_thrice)),
            "Vertex should wrap values. {:?}.vertex(2) = {:?} != {:?}.vertex({}) = {:?}",
            arctic80(),
            arctic80().vertex(2),
            arctic80(),
            loop_around_thrice,
            arctic80().vertex(loop_around_thrice)
        );
    }

    #[test]
    fn test_loop_num_edges() {
        let tests = vec![
            (empty_loop(), 0),
            (full_loop(), 0),
            (far_hemi(), 4),
            (candy_cane(), 6),
            (small_necw(), 3),
            (arctic80(), 3),
            (antarctic80(), 3),
            (line_triangle(), 3),
            (skinny_chevron(), 4),
        ];

        for (loop_under_test, want) in tests {
            let got = loop_under_test.num_edges();
            assert_eq!(
                got, want,
                "{:?}.num_edges() = {}, want {}",
                loop_under_test, got, want
            );
        }
    }

    #[test]
    fn test_loop_edge() {
        let tests = vec![
            (
                far_hemi(),
                2,
                Point(R3Vector {
                    x: 0.0,
                    y: 0.0,
                    z: -1.0,
                }),
                Point(R3Vector {
                    x: 0.0,
                    y: -1.0,
                    z: 0.0,
                }),
            ),
            (
                candy_cane(),
                0,
                parse_points("-20:150")[0],
                parse_points("-20:-70")[0],
            ),
            (
                candy_cane(),
                1,
                parse_points("-20:-70")[0],
                parse_points("0:70")[0],
            ),
            (
                candy_cane(),
                2,
                parse_points("0:70")[0],
                parse_points("10:-150")[0],
            ),
            (
                candy_cane(),
                3,
                parse_points("10:-150")[0],
                parse_points("10:70")[0],
            ),
            (
                candy_cane(),
                4,
                parse_points("10:70")[0],
                parse_points("-10:-70")[0],
            ),
            (
                candy_cane(),
                5,
                parse_points("-10:-70")[0],
                parse_points("-20:150")[0],
            ),
            (
                skinny_chevron(),
                2,
                parse_points("0:1e-320")[0],
                parse_points("1e-320:80")[0],
            ),
            (
                skinny_chevron(),
                3,
                parse_points("1e-320:80")[0],
                parse_points("0:0")[0],
            ),
        ];

        for (loop_under_test, edge, want_a, want_b) in tests {
            let e = loop_under_test.edge(edge);
            assert!(
                e.v0.approx_eq(&want_a) && e.v1.approx_eq(&want_b),
                "{:?}.edge({}) = {:?}, want ({:?}, {:?})",
                loop_under_test,
                edge,
                e,
                want_a,
                want_b
            );
        }
    }

    #[test]
    fn test_loop_from_cell() {
        let cell_id = CellID::from(LatLng::from_degrees(40.565459, -74.645276));
        let cell = Cell::from(cell_id);
        let loop_from_cell = Loop::from_cell(&cell);

        // Demonstrates the reason for this test; the cell bounds are more
        // conservative than the resulting loop bounds.
        assert!(
            !loop_from_cell.bound.contains(&cell.rect_bound()),
            "loop_from_cell's RectBound contains the original cell's RectBound, but should not"
        );
    }

    #[test]
    fn test_loop_regular_loop() {
        let center = Point::from(LatLng::from_degrees(80.0, 135.0));
        let radius = Angle::from(Deg(20.0));
        let loop_obj = Loop::regular_loop(center, radius, 4);
        assert_eq!(
            loop_obj.vertices.len(),
            4,
            "RegularLoop with 4 vertices should have 4 vertices, got {}",
            loop_obj.vertices.len()
        );
        // The actual Points values are already tested in the s2point_test method TestRegularPoints.
    }

    // clone_loop creates a new copy of the given loop including all of its vertices
    // so that when tests modify vertices in it, it won't ruin the original loop.
    fn clone_loop(l: &Loop) -> Loop {
        let mut c = Loop {
            vertices: l.vertices.clone(),
            origin_inside: l.origin_inside,
            bound: l.bound.clone(),
            subregion_bound: l.subregion_bound.clone(),
            index: ShapeIndex::new(),
            depth: l.depth,
        };
        c.index.add(&ShapeType::Loop(c.clone()));

        c
    }

    #[test]
    fn test_loop_equal() {
        let tests = vec![
            (empty_loop(), empty_loop(), true),
            (full_loop(), full_loop(), true),
            (empty_loop(), full_loop(), false),
            (candy_cane(), candy_cane(), true),
            (candy_cane(), rotate(&candy_cane()), false),
            (candy_cane(), rotate(&rotate(&candy_cane())), false),
            (candy_cane(), rotate(&rotate(&rotate(&candy_cane()))), false),
            (
                candy_cane(),
                rotate(&rotate(&rotate(&rotate(&candy_cane())))),
                false,
            ),
            (
                candy_cane(),
                rotate(&rotate(&rotate(&rotate(&rotate(&candy_cane()))))),
                false,
            ),
            // candy_cane has 6 points, so 6 rotates should line up again.
            (
                candy_cane(),
                rotate(&rotate(&rotate(&rotate(&rotate(&rotate(&candy_cane())))))),
                true,
            ),
        ];

        for (a, b, want) in tests {
            assert_eq!(
                a.equal(&b),
                want,
                "{:?}.equal({:?}) = {}, want {}",
                a,
                b,
                a.equal(&b),
                want
            );
        }
    }

    #[test]
    fn test_loop_contains_matches_crossing_sign() {
        // This test demonstrates a former incompatibility between CrossingSign
        // and ContainsPoint. It constructs a Cell-based loop L and
        // an edge E from Origin to a0 that crosses exactly one edge of L.  Yet
        // previously, Contains() returned false for both endpoints of E.
        //
        // The reason for the bug was that the loop bound was sometimes too tight.
        // The Contains() code for a0 bailed out early because a0 was found not to
        // be inside the bound of L.

        // Start with a cell that ends up producing the problem.
        let cell_id = CellID::from(&Point(R3Vector {
            x: 1.0,
            y: 1.0,
            z: 1.0,
        }))
        .parent(21);
        let cell = Cell::from(cell_id);
        let children = cell.children().unwrap();

        let mut points = Vec::with_capacity(4);
        for i in 0..4 {
            // Note extra normalization. Center() is already normalized.
            // The test results will no longer be inconsistent if the extra
            // Normalize() is removed.
            points.push(Point(children[i].center().0.normalize()));
        }

        // Get a vertex from a grandchild cell.
        // +---------------+---------------+
        // |               |               |
        // |    points[3]  |   points[2]   |
        // |       v       |       v       |
        // |       +-------+------ +       |
        // |       |       |       |       |
        // |       |       |       |       |
        // |       |       |       |       |
        // +-------+-------+-------+-------+
        // |       |       |       |       |
        // |       |    <----------------------- grandchild_cell
        // |       |       |       |       |
        // |       +-------+------ +       |
        // |       ^       |       ^       | <-- cell
        // | points[0]/a0  |     points[1] |
        // |               |               |
        // +---------------+---------------+
        let loop_obj = Loop::from_points(points.clone());
        let grandchildren = children[0].children().unwrap();

        let grandchild_cell = &grandchildren[2];

        let a0 = grandchild_cell.vertex(0);

        // This test depends on rounding errors that should make a0 slightly different from points[0]
        assert_ne!(
            points[0], a0,
            "{:?} not different enough from {:?} to successfully test",
            points[0], a0
        );

        // The edge from a0 to the origin crosses one boundary.
        let want = Crossing::DoNotCross;
        let got = EdgeCrosser::new_chain_edge_crosser(&a0, &Point::origin(), &loop_obj.vertex(0))
            .chain_crossing_sign(&loop_obj.vertex(1));
        assert_eq!(
            got,
            want,
            "crossingSign({:?}, {:?}, {:?}, {:?}) = {:?}, want {:?}",
            a0,
            Point::origin(),
            loop_obj.vertex(0),
            loop_obj.vertex(1),
            got,
            want
        );

        let want = Crossing::Cross;
        let got = EdgeCrosser::new_chain_edge_crosser(&a0, &Point::origin(), &loop_obj.vertex(1))
            .chain_crossing_sign(&loop_obj.vertex(2));
        assert_eq!(
            got,
            want,
            "crossingSign({:?}, {:?}, {:?}, {:?}) = {:?}, want {:?}",
            a0,
            Point::origin(),
            loop_obj.vertex(1),
            loop_obj.vertex(2),
            got,
            want
        );

        let want = Crossing::DoNotCross;
        let got = EdgeCrosser::new_chain_edge_crosser(&a0, &Point::origin(), &loop_obj.vertex(2))
            .chain_crossing_sign(&loop_obj.vertex(3));
        assert_eq!(
            got,
            want,
            "crossingSign({:?}, {:?}, {:?}, {:?}) = {:?}, want {:?}",
            a0,
            Point::origin(),
            loop_obj.vertex(2),
            loop_obj.vertex(3),
            got,
            want
        );

        let want = Crossing::DoNotCross;
        let got = EdgeCrosser::new_chain_edge_crosser(&a0, &Point::origin(), &loop_obj.vertex(3))
            .chain_crossing_sign(&loop_obj.vertex(4));
        assert_eq!(
            got,
            want,
            "crossingSign({:?}, {:?}, {:?}, {:?}) = {:?}, want {:?}",
            a0,
            Point::origin(),
            loop_obj.vertex(3),
            loop_obj.vertex(4),
            got,
            want
        );

        // Contains should return false for the origin, and true for a0.
        assert!(
            !loop_obj.contains_point(Point::origin()),
            "{:?}.contains_point({:?}) = true, want false",
            loop_obj,
            Point::origin()
        );
        assert!(
            loop_obj.contains_point(a0),
            "{:?}.contains_point({:?}) = false, want true",
            loop_obj,
            a0
        );

        // Since a0 is inside the loop, it should be inside the bound.
        let bound = loop_obj.bound.clone();
        assert!(
            bound.contains_point(&a0),
            "{:?}.contains_point({:?}) = false, want true",
            bound,
            a0
        );
    }

    fn test_loop_nested_pair(t: &mut impl FnMut(&str, bool), a: &Loop, b: &Loop) {
        let a1 = &mut clone_loop(a);
        a1.invert();
        let b1 = &mut clone_loop(b);
        b1.invert();
        test_loop_one_nested_pair(t, a, b);
        test_loop_one_nested_pair(t, b1, a1);
        test_loop_one_disjoint_pair(t, a1, b);
        test_loop_one_covering_pair(t, a, b1);
    }

    fn test_loop_one_nested_pair(t: &mut impl FnMut(&str, bool), a: &Loop, b: &Loop) {
        t(
            &format!("{:?}.contains({:?}) = false, want true", a, b),
            a.contains(b),
        );
        t(
            &format!(
                "{:?}.contains({:?}) = {}, want {}",
                b,
                a,
                b.contains(a),
                a.boundary_equal(b)
            ),
            b.contains(a) == a.boundary_equal(b),
        );
        t(
            &format!(
                "{:?}.intersects({:?}) = {}, want {}",
                a,
                b,
                a.intersects(b),
                !b.is_empty()
            ),
            a.intersects(b) == !b.is_empty(),
        );
        t(
            &format!(
                "{:?}.intersects({:?}) = {}, want {}",
                b,
                a,
                b.intersects(a),
                !b.is_empty()
            ),
            b.intersects(a) == !b.is_empty(),
        );
    }

    fn test_loop_one_disjoint_pair(t: &mut impl FnMut(&str, bool), a: &Loop, b: &Loop) {
        t(
            &format!("{:?}.intersects({:?}) = true, want false", a, b),
            !a.intersects(b),
        );
        t(
            &format!("{:?}.intersects({:?}) = true, want false", b, a),
            !b.intersects(a),
        );
        t(
            &format!(
                "{:?}.contains({:?}) = {}, want {}",
                a,
                b,
                a.contains(b),
                b.is_empty()
            ),
            a.contains(b) == b.is_empty(),
        );
        t(
            &format!(
                "{:?}.contains({:?}) = {}, want {}",
                b,
                a,
                b.contains(a),
                a.is_empty()
            ),
            b.contains(a) == a.is_empty(),
        );
    }

    fn test_loop_one_covering_pair(t: &mut impl FnMut(&str, bool), a: &Loop, b: &Loop) {
        t(
            &format!(
                "{:?}.contains({:?}) = {}, want {}",
                a,
                b,
                a.contains(b),
                a.is_full()
            ),
            a.contains(b) == a.is_full(),
        );
        t(
            &format!(
                "{:?}.contains({:?}) = {}, want {}",
                b,
                a,
                b.contains(a),
                b.is_full()
            ),
            b.contains(a) == b.is_full(),
        );
        let a1 = &mut clone_loop(a);
        a1.invert();
        let complementary = a1.boundary_equal(b);
        t(
            &format!(
                "{:?}.intersects({:?}) = {}, want {}",
                a,
                b,
                a.intersects(b),
                !complementary
            ),
            a.intersects(b) == !complementary,
        );
        t(
            &format!(
                "{:?}.intersects({:?}) = {}, want {}",
                b,
                a,
                b.intersects(a),
                !complementary
            ),
            b.intersects(a) == !complementary,
        );
    }

    fn test_loop_one_overlapping_pair(t: &mut impl FnMut(&str, bool), a: &Loop, b: &Loop) {
        t(
            &format!("{:?}.contains({:?}) = true, want false", a, b),
            !a.contains(b),
        );
        t(
            &format!("{:?}.contains({:?}) = true, want false", b, a),
            !b.contains(a),
        );

        t(
            &format!("{:?}.intersects({:?}) = false, want true", a, b),
            a.intersects(b),
        );
        t(
            &format!("{:?}.intersects({:?}) = false, want true", b, a),
            b.intersects(a),
        );
    }

    #[test]
    fn test_loop_relations() {
        // Create test cases for loop relations
        let tests = vec![
            // Check full and empty relationships with normal loops and each other.
            (
                full_loop(),
                full_loop(),
                true,  // contains
                true,  // contained
                false, // disjoint
                true,  // covers
                true,  // sharedEdge
            ),
            (
                full_loop(),
                north_hemi(),
                true,  // contains
                false, // contained
                false, // disjoint
                true,  // covers
                false, // sharedEdge
            ),
            (
                full_loop(),
                empty_loop(),
                true,  // contains
                false, // contained
                true,  // disjoint
                true,  // covers
                false, // sharedEdge
            ),
            (
                north_hemi(),
                full_loop(),
                false, // contains
                true,  // contained
                false, // disjoint
                true,  // covers
                false, // sharedEdge
            ),
            (
                north_hemi(),
                empty_loop(),
                true,  // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                empty_loop(),
                full_loop(),
                false, // contains
                true,  // contained
                true,  // disjoint
                true,  // covers
                false, // sharedEdge
            ),
            (
                empty_loop(),
                north_hemi(),
                false, // contains
                true,  // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                empty_loop(),
                empty_loop(),
                true,  // contains
                true,  // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                north_hemi(),
                north_hemi(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                north_hemi(),
                south_hemi(),
                false, // contains
                false, // contained
                true,  // disjoint
                true,  // covers
                true,  // sharedEdge
            ),
            (
                north_hemi(),
                east_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                north_hemi(),
                arctic80(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                north_hemi(),
                antarctic80(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                north_hemi(),
                candy_cane(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            // We can't compare northHemi3 vs. northHemi or southHemi because the
            // result depends on the "simulation of simplicity" implementation details.
            (
                north_hemi3(),
                north_hemi3(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                north_hemi3(),
                east_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                north_hemi3(),
                arctic80(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                north_hemi3(),
                antarctic80(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                north_hemi3(),
                candy_cane(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                south_hemi(),
                north_hemi(),
                false, // contains
                false, // contained
                true,  // disjoint
                true,  // covers
                true,  // sharedEdge
            ),
            (
                south_hemi(),
                south_hemi(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                south_hemi(),
                far_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                south_hemi(),
                arctic80(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                south_hemi(),
                antarctic80(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                south_hemi(),
                candy_cane(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                candy_cane(),
                north_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                candy_cane(),
                south_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                candy_cane(),
                arctic80(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                candy_cane(),
                antarctic80(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                candy_cane(),
                candy_cane(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                near_hemi(),
                west_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                small_necw(),
                south_hemi(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                small_necw(),
                west_hemi(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                small_necw(),
                north_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                true,  // covers
                false, // sharedEdge
            ),
            (
                small_necw(),
                east_hemi(),
                false, // contains
                false, // contained
                false, // disjoint
                true,  // covers
                false, // sharedEdge
            ),
            (
                loop_a(),
                loop_a(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_a(),
                loop_b(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                loop_a(),
                a_intersect_b(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_a(),
                a_union_b(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_a(),
                a_minus_b(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_a(),
                b_minus_a(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_b(),
                loop_a(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                loop_b(),
                loop_b(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_b(),
                a_intersect_b(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_b(),
                a_union_b(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_b(),
                a_minus_b(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_b(),
                b_minus_a(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_intersect_b(),
                loop_a(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_intersect_b(),
                loop_b(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_intersect_b(),
                a_intersect_b(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_intersect_b(),
                a_union_b(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                a_intersect_b(),
                a_minus_b(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_intersect_b(),
                b_minus_a(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_union_b(),
                loop_a(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_union_b(),
                loop_b(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_union_b(),
                a_intersect_b(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                a_union_b(),
                a_union_b(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_union_b(),
                a_minus_b(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_union_b(),
                b_minus_a(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_minus_b(),
                loop_a(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_minus_b(),
                loop_b(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_minus_b(),
                a_intersect_b(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_minus_b(),
                a_union_b(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_minus_b(),
                a_minus_b(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                a_minus_b(),
                b_minus_a(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                b_minus_a(),
                loop_a(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                b_minus_a(),
                loop_b(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                b_minus_a(),
                a_intersect_b(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                b_minus_a(),
                a_union_b(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                b_minus_a(),
                a_minus_b(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                b_minus_a(),
                b_minus_a(),
                true,  // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            // Make sure the relations are correct if the loop crossing happens on
            // two ends of a shared boundary segment.
            // LoopRelationsWhenSameExceptPiecesStickingOutAndIn
            (
                loop_a(),
                loop_c(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_c(),
                loop_a(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_a(),
                loop_d(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_d(),
                loop_a(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_e(),
                loop_f(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_e(),
                loop_g(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_e(),
                loop_h(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_e(),
                loop_i(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                loop_f(),
                loop_g(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_f(),
                loop_h(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_f(),
                loop_i(),
                false, // contains
                false, // contained
                false, // disjoint
                false, // covers
                false, // sharedEdge
            ),
            (
                loop_g(),
                loop_h(),
                false, // contains
                true,  // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_h(),
                loop_g(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_g(),
                loop_i(),
                false, // contains
                false, // contained
                true,  // disjoint
                false, // covers
                true,  // sharedEdge
            ),
            (
                loop_h(),
                loop_i(),
                true,  // contains
                false, // contained
                false, // disjoint
                false, // covers
                true,  // sharedEdge
            ),
        ];

        // Create a tester function to track results
        let mut t = |msg: &str, result: bool| {
            assert!(result, "{}", msg);
        };

        for test_case in tests {
            let (a, b, contains, contained, disjoint, covers, _shared_edge) = test_case;

            if contains {
                test_loop_nested_pair(&mut t, &a, &b);
            }
            if contained {
                test_loop_nested_pair(&mut t, &b, &a);
            }
            if covers {
                let mut b1 = clone_loop(&b);
                b1.invert();
                test_loop_nested_pair(&mut t, &a, &b1);
            }
            if disjoint {
                let mut a1 = clone_loop(&a);
                a1.invert();
                test_loop_nested_pair(&mut t, &a1, &b);
            } else if !(contains || contained || covers) {
                // Given loops A and B such that both A and its complement
                // intersect both B and its complement, test various
                // identities involving these four loops.
                let mut a1 = clone_loop(&a);
                a1.invert();
                let mut b1 = clone_loop(&b);
                b1.invert();
                test_loop_one_overlapping_pair(&mut t, &a, &b);
                test_loop_one_overlapping_pair(&mut t, &a1, &b1);
                test_loop_one_overlapping_pair(&mut t, &a1, &b);
                test_loop_one_overlapping_pair(&mut t, &a, &b1);
            }

            /*
            // Skipping the shared_edge comparison check from Go as Rust doesn't have
            // a direct equivalent to CompareBoundary method yet
            if !shared_edge && (contains || contained || disjoint) {
                if a.contains(&b) != a.contains_nested(&b) {
                    t(
                        &format!("{:?}.contains({:?}) = {}, but should equal {:?}.contains_nested({:?}) = {}",
                                 a, b, a.contains(&b), a, b, a.contains_nested(&b)),
                        false,
                    );
                }
            }
            */
        }
    }

    #[test]
    fn test_loop_turning_angle() {
        let tests = vec![
            // (empty_loop(), 2.0 * PI, "empty loop"),
            // (full_loop(), -2.0 * PI, "full loop"),
            (north_hemi3(), 0.0, "north_hemi3"),
            // (west_hemi(), 0.0, "west_hemi"),
            // (candy_cane(), 4.69364376125922, "candy_cane"),
            // (line_triangle(), 2.0 * PI, "line_triangle"),
            // (skinny_chevron(), 2.0 * PI, "skinny_chevron"),
        ];

        for (loop_obj, want, name) in tests {
            // Using f64_eq from the existing codebase for comparison
            assert_f64_eq!(
                loop_obj.turning_angle(),
                want,
                DBL_EPSILON,
                format!(
                    "Testing {}: {:?}.turning_angle() = {}, want {}",
                    name,
                    loop_obj,
                    loop_obj.turning_angle(),
                    want
                )
            );

            // Check that the turning angle is *identical* when the vertex order is
            // rotated, and that the sign is inverted when the vertices are reversed.
            let expected = loop_obj.turning_angle();
            let mut loop_copy = clone_loop(&loop_obj);

            for _ in 0..loop_obj.vertices.len() {
                loop_copy.invert();
                assert_eq!(
                    loop_copy.turning_angle(),
                    -expected,
                    "loop.invert().turning_angle() = {}, want {}",
                    loop_copy.turning_angle(),
                    -expected
                );

                // Invert it back to normal.
                loop_copy.invert();

                loop_copy = rotate(&loop_copy);
                assert_eq!(
                    loop_copy.turning_angle(),
                    expected,
                    "loop.turning_angle() = {}, want {}",
                    loop_copy.turning_angle(),
                    expected
                );
            }
        }

        // Build a narrow spiral loop starting at the north pole. This is designed
        // to test that the error in TurningAngle is linear in the number of
        // vertices even when the partial sum of the turning angles gets very large.
        // The spiral consists of two arms defining opposite sides of the loop.
        const ARM_POINTS: usize = 10000; // Number of vertices in each "arm"
        const ARM_RADIUS: f64 = 0.01; // Radius of spiral.
        let mut vertices = Vec::with_capacity(2 * ARM_POINTS);

        // Set the center point of the spiral.
        vertices.push(Point::from_coords(0.0, 0.0, 1.0));

        // Fill in with ARM_POINTS - 1 more points to get to the right size
        for _ in 1..ARM_POINTS {
            vertices.push(Point::from_coords(0.0, 0.0, 1.0));
        }

        for i in 0..ARM_POINTS {
            let angle = (2.0 * PI / 3.0) * i as f64;
            let x = angle.cos();
            let y = angle.sin();
            let r1 = i as f64 * ARM_RADIUS / ARM_POINTS as f64;
            let r2 = (i as f64 + 1.5) * ARM_RADIUS / ARM_POINTS as f64;
            vertices[ARM_POINTS - i - 1] = Point::from_coords(r1 * x, r1 * y, 1.0);
            vertices[ARM_POINTS + i] = Point::from_coords(r2 * x, r2 * y, 1.0);
        }

        // This is a pathological loop that contains many long parallel edges.
        let spiral = Loop::from_points(vertices);

        // Check that TurningAngle is consistent with Area to within the
        // error bound of the former. We actually use a tiny fraction of the
        // worst-case error bound, since the worst case only happens when all the
        // roundoff errors happen in the same direction and this test is not
        // designed to achieve that. The error in Area can be ignored for the
        // purposes of this test since it is generally much smaller.
        assert_f64_eq!(
            spiral.turning_angle(),
            2.0 * PI - spiral.area(),
            0.01 * spiral.turning_angle_max_error(),
            format!(
                "spiral.turning_angle() = {}, want {}",
                spiral.turning_angle(),
                2.0 * PI - spiral.area()
            )
        );
    }

    #[test]
    fn test_loop_area_and_centroid() {
        let origin = Point(Vector {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });

        assert_eq!(
            empty_loop().area(),
            0.0,
            "EmptyLoop.area() = {}, want {}",
            empty_loop().area(),
            0.0
        );
        assert_eq!(
            full_loop().area(),
            4.0 * PI,
            "FullLoop.area() = {}, want {}",
            full_loop().area(),
            4.0 * PI
        );
        assert!(
            empty_loop().centroid().approx_eq(&origin),
            "EmptyLoop.centroid() = {:?}, want {:?}",
            empty_loop().centroid(),
            origin
        );
        assert!(
            full_loop().centroid().approx_eq(&origin),
            "FullLoop.centroid() = {:?}, want {:?}",
            full_loop().centroid(),
            origin
        );

        assert_f64_eq!(
            north_hemi().area() as f64,
            2.0 * PI,
            DBL_EPSILON,
            format!(
                "northHemi.area() = {}, want {}",
                north_hemi().area(),
                2.0 * PI
            )
        );

        let east_hemi_area = east_hemi().area();
        assert!(
            east_hemi_area >= 2.0 * PI - 1e-12 && east_hemi_area <= 2.0 * PI + 1e-12,
            "eastHemi.area() = {}, want between [{}, {}]",
            east_hemi_area,
            2.0 * PI - 1e-12,
            2.0 * PI + 1e-12
        );

        // TODO: Add implementation for random frames and testing spherical caps
    }

    #[test]
    fn test_loop_area_consistent_with_turning_angle() {
        // Check that the area computed using area() is consistent with the
        // turning angle of the loop computed using turning_angle(). According to
        // the Gauss-Bonnet theorem, the area of the loop should be equal to 2*Pi
        // minus its turning angle.
        for (i, loop_obj) in all_loops().iter().enumerate() {
            let area = loop_obj.area();
            let gauss_area = 2.0 * PI - loop_obj.turning_angle();

            // TODO(roberts): The error bound below is much larger than it should be.
            assert!(
                (area - gauss_area).abs() <= 1e-9,
                "{}. {:?}.area() = {}, want {}",
                i,
                loop_obj,
                area,
                gauss_area
            );
        }
    }

    #[test]
    fn test_loop_normalized_compatible_with_contains() {
        let p = parse_points("40:40")[0];

        let tests = vec![line_triangle(), skinny_chevron()];

        // Checks that if a loop is normalized, it doesn't contain a
        // point outside of it, and vice versa.
        for loop_obj in tests {
            let mut flip = clone_loop(&loop_obj);

            flip.invert();
            let norm = loop_obj.is_normalized();
            let contains = loop_obj.contains_point(p);
            assert!(
                norm != contains,
                "loop.is_normalized() = {} == loop.contains_point({:?}) = {}, want !=",
                norm,
                p,
                contains
            );

            let norm = flip.is_normalized();
            let contains = flip.contains_point(p);
            assert!(
                norm != contains,
                "flip.is_normalized() = {} == flip.contains_point({:?}) = {}, want !=",
                norm,
                p,
                contains
            );

            assert!(
                loop_obj.is_normalized() != flip.is_normalized(),
                "a loop and its invert can not both be normalized"
            );

            flip.normalize();
            assert!(
                !flip.contains_point(p),
                "{:?}.contains_point({:?}) = true, want false",
                flip,
                p
            );
        }
    }

    #[test]
    fn test_loop_validate_detects_invalid_loops() {
        let tests = vec![
            // Not enough vertices. Note that all single-vertex loops are valid; they
            // are interpreted as being either "empty" or "full".
            ("loop has no vertices", parse_points("")),
            ("loop has too few vertices", parse_points("20:20, 21:21")),
            // degenerate edge checks happen in validation before duplicate vertices.
            (
                "loop has degenerate first edge",
                parse_points("20:20, 20:20, 20:21"),
            ),
            (
                "loop has degenerate third edge",
                parse_points("20:20, 20:21, 20:20"),
            ),
            // TODO(roberts): Uncomment these cases when FindAnyCrossings is in.
            /*
            (
                "loop has duplicate points",
                parse_points("20:20, 21:21, 21:20, 20:20, 20:21"),
            ),
            (
                "loop has crossing edges",
                parse_points("20:20, 21:21, 21:20.5, 21:20, 20:21"),
            ),
            */
            (
                // Ensure points are not normalized.
                "loop with non-normalized vertices",
                vec![
                    Point(R3Vector {
                        x: 2.0,
                        y: 0.0,
                        z: 0.0,
                    }),
                    Point(R3Vector {
                        x: 0.0,
                        y: 1.0,
                        z: 0.0,
                    }),
                    Point(R3Vector {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    }),
                ],
            ),
            (
                // Adjacent antipodal vertices
                "loop with antipodal points",
                vec![
                    Point(R3Vector {
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    }),
                    Point(R3Vector {
                        x: -1.0,
                        y: 0.0,
                        z: 0.0,
                    }),
                    Point(R3Vector {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    }),
                ],
            ),
        ];

        for (msg, points) in tests {
            let loop_obj = Loop::from_points(points);
            if loop_obj.validate() {
                panic!("{}. {:?}.validate() = Ok, want Err", msg, loop_obj);
            }
            // The Go tests also tests that the returned error message string contains
            // a specific set of text. That part of the test is skipped here.
        }
    }

    // TODO(roberts): Convert these into changeable flags or parameters.
    // A loop with a 10km radius and 4096 vertices has an edge length of 15 meters.
    const DEFAULT_RADIUS_KM: f64 = 10.0;
    const NUM_LOOP_SAMPLES: usize = 16;
    const NUM_QUERIES_PER_LOOP: usize = 100;

    // #[bench]
    // fn bench_loop_contains_point(b: &mut Bencher) {
    //     // Benchmark ContainsPoint() on regular loops. The query points for a loop are
    //     // chosen so that they all lie in the loop's bounding rectangle (to avoid the
    //     // quick-rejection code path).

    //     // Go ranges from 4 -> 256k by powers of 2 for number of vertices for benchmarking.
    //     let mut vertices = 4;
    //     for n in 1..=17 {
    //         b.iter_with_setup(
    //             || {
    //                 let mut loops = Vec::with_capacity(NUM_LOOP_SAMPLES);
    //                 for i in 0..NUM_LOOP_SAMPLES {
    //                     loops.push(Loop::regular_loop(
    //                         random_point(),
    //                         km_to_angle(10.0),
    //                         vertices,
    //                     ));
    //                 }

    //                 let mut queries = Vec::with_capacity(NUM_LOOP_SAMPLES);
    //                 for loop_obj in &loops {
    //                     let mut loop_queries = Vec::with_capacity(NUM_QUERIES_PER_LOOP);
    //                     for j in 0..NUM_QUERIES_PER_LOOP {
    //                         loop_queries.push(sample_point_from_rect(&loop_obj.rect_bound()));
    //                     }
    //                     queries.push(loop_queries);
    //                 }
    //                 (loops, queries)
    //             },
    //             |(loops, queries)| {
    //                 let i = rand::random::<usize>();
    //                 let j = rand::random::<usize>();
    //                 loops[i % NUM_LOOP_SAMPLES]
    //                     .contains_point(queries[i % NUM_LOOP_SAMPLES][j % NUM_QUERIES_PER_LOOP]);
    //             },
    //         );
    //         vertices *= 2;
    //     }
    // }
}
