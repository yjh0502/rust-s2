#[allow(dead_code)]
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

// This file contains a collection of methods for:
//
//   (1) Robustly clipping geodesic edges to the faces of the S2 biunit cube
//       (see s2stuv), and
//
//   (2) Robustly clipping 2D edges against 2D rectangles.
//
// These functions can be used to efficiently find the set of CellIDs that
// are intersected by a geodesic edge (e.g., see crossing_edge_query).

use consts::*;
use point::Point;
use r2;
use r3;
use r1;
use s2::stuv;

// EDGE_CLIP_ERROR_UV_COORD is the maximum error in a u- or v-coordinate compared to the exact result, assuming that the points A and B are in
// the rectangle [-1,1]x[1,1] or slightly outside it (by 1e-10 or less).
#[allow(dead_code)]
const EDGE_CLIP_ERROR_UV_COORD: f64 = 2.25 * DBL_EPSILON;

// EDGE_CLIP_ERORR_UV_DIST is the maximum distance from a clipped point to
// the corresponding exact result. It is equal to the error in a single
// coordinate because at most one coordinate is subject to error.
#[allow(dead_code)]
const EDGE_CLIP_ERROR_UV_DIST: f64 = 2.25 * DBL_EPSILON;

// faceClipErrorRadians is the maximum angle between a returned vertex
// and the nearest point on the exact edge AB. It is equal to the
// maximum directional error in PointCross, plus the error when
// projecting points onto a cube face
const FACE_CLIP_ERROR_RADIANS: f64 = 3.0 * DBL_EPSILON;

// faceClipErrorDist is the same angle expressed as a maximum distance
// in (u,v)-space. In other words, a returned vertex is at most this far
// from the exact edge AB projected into (u,v)-space.
#[allow(dead_code)]
const FACE_CLIP_ERROR_UV_DIST: f64 = 9.0 * DBL_EPSILON;

// faceClipErrorUVCoord is the maximum angle between a returned vertex
// and the nearest point on the exact edge AB expressed as the maximum error
// in an individual u- or v-coordinate. In other words, for each
// returned vertex there is a point on the exact edge AB whose u- and
// v-coordinates differ from the vertex by at most this amount
const FACE_CLIP_ERROR_UV_COORD: f64 = 9.0 * (1.0 / std::f64::consts::SQRT_2) * DBL_EPSILON;

// intersectsRectErrorUVDist is the maximum error when computing if a point
// intersects with a given Rect. If some point of AB is inside the
// rectangle by at least this distance, the result is guaranteed to be true;
// if all points of AB are outside the rectangle by at least this distance,
// the result is guaranteed to be false. This bound assumes that rect is
// a subset of the rectangle [-1,1]x[-1,1] or extends slightly outside it
// (e.g., by 1e-10 or less)
#[allow(dead_code)]
const INTERSECT_RECT_ERROR_UV_DIST: f64 = 3.0 * std::f64::consts::SQRT_2 * DBL_EPSILON;

// clip_to_face returns the (u,v) coordinates for the portion of the edge AB that
// intersects the given face, or false if the edge AB does not intersect.
// This method guarantees that the clipped vertices lie within the [-1,1]x[-1,1]
// cube face rectangle and are within faceClipErrorUVDist of the line AB, but
// the results may differ from those produced by FaceSegments.
pub fn clip_to_face(a: Point, b: Point, face: u8) -> (r2::point::Point, r2::point::Point, bool) {
    clip_to_padded_face(a, b, face, 0.0)
}

// clip_to_padded_face returns the (u,v) coordinates for the portion of the edge AB that
// intersects the given face, but rather than clipping to the square [-1,1]x[-1,1]
// in (u,v) space, this method clips to [-R,R]x[-R,R] where R=(1+padding).
// Padding must be non-negative.
pub fn clip_to_padded_face(a: Point, b: Point, f: u8, padding: f64) -> (r2::point::Point, r2::point::Point, bool) {
    // Fast path: both endpoints are on the given face
    if stuv::face(&a.0) == f && stuv::face(&b.0) == f {
        let (au, av) = stuv::valid_face_xyz_to_uv(f, &a.0);
        let (bu, bv) = stuv::valid_face_xyz_to_uv(f, &a.0);
        return (r2::point::Point{x: au, y: av}, r2::point::Point{x: bu, y: bv}, true)
    }

    // Convert everything into the (u,v,w) coordinates of the given face. Note
	// that the cross product *must* be computed in the original (x,y,z)
	// coordinate system because PointCross (unlike the mathematical cross
	// product) can produce different results in different coordinate systems
	// when one argument is a linear multiple of the other, due to the use of
	// symbolic perturbations.
    let mut norm_uvw: PointUVW = stuv::face_xyz_to_uvw(f, &a.cross(&b));
    let a_uvw: PointUVW = stuv::face_xyz_to_uvw(f, &a);
    let b_uvw: PointUVW = stuv::face_xyz_to_uvw(f, &b);

	// Padding is handled by scaling the u- and v-components of the normal.
	// Letting R=1+padding, this means that when we compute the dot product of
	// the normal with a cube face vertex (such as (-1,-1,1)), we will actually
	// compute the dot product with the scaled vertex (-R,-R,1). This allows
	// methods such as intersectsFace, exitAxis, etc, to handle padding
	// with no further modifications.
    let scale_uv = 1.0 + padding;
    let scaled_n: PointUVW = Point(r3::vector::Vector{x: scale_uv * norm_uvw.0.x, y: scale_uv * norm_uvw.0.y, z: norm_uvw.0.z});
    if !scaled_n.intersects_face(){ 
        return (r2::point::Point{x: 0.0, y:0.0}, r2::point::Point{x: 0.0, y:0.0}, false)
     }
    let ldexp = |x: f64, exp: f64| -> f64 {
        x * (exp.exp2() as f64) 
    };

    // TODO: This is a workaround for extremely small vectors where some
	// loss of precision can occur in Normalize causing underflow. When PointCross
	// is updated to work around this, this can be removed.
    if norm_uvw.0.x.abs().max(norm_uvw.0.y.abs().max(norm_uvw.0.z.abs())) < ldexp(1.0, -511.0) {
        norm_uvw = norm_uvw * ldexp(1.0, 563.0)
    }
    norm_uvw = norm_uvw.normalize();
    let a_tan = norm_uvw.cross(&a_uvw);
    let b_tan = norm_uvw.cross(&b_uvw);

    // As described in clipDestination, if the sum of the scores from clipping the two
	// endpoints is 3 or more, then the segment does not intersect this face
    let (a_uv, a_score) = clip_destination(b_uvw, a_uvw, scaled_n * -1_f64, b_tan, a_tan, scale_uv);
    let (b_uv, b_score) = clip_destination(a_uvw, b_uvw, scaled_n * -1_f64, a_tan, b_tan, scale_uv);

    (a_uv, b_uv, a_score+b_score < 3)
}
 
// clip_edge returns the portion of the edge defined by AB that is contained by the
// given rectangle. If there is no intersection, false is returned and aClip and bClip
// are undefined.
pub fn clip_edge(a: r2::point::Point, b: r2::point::Point, clip: r2::rect::Rect) -> (r2::point::Point, r2::point::Point, bool) {
    // Compute the bounding rectangle of AB, clip it, and then extract the new
	// endpoints from the clipped bound
    let bound = r2::rect::Rect::from_points(&[a,b]);
    let (bound, intersects) = clip_edge_bound(a, b, clip, bound); 
    if !intersects {
        return (
            r2::point::Point{x: 0.0, y: 0.0},
            r2::point::Point{x: 0.0, y: 0.0},
            intersects
        )
    }
    let ai = if a.x > b.x {
        1
    } else {
        0
    };

    let aj = if a.y > b.y {
        1
    } else {
        0
    };
    (bound.vertex_ij(ai, aj), bound.vertex_ij(1-ai, 1-aj), true)
}

// The three functions below (sumEqual, intersectsFace, intersectsOppositeEdges)
// all compare a sum (u + v) to a third value w. They are implemented in such a
// way that they produce an exact result even though all calculations are done
// with ordinary floating-point operations. Here are the principles on which these
// functions are based:
//
// A. If u + v < w in floating-point, then u + v < w in exact arithmetic.
//
// B. If u + v < w in exact arithmetic, then at least one of the following
//    expressions is true in floating-point:
//       u + v < w
//       u < w - v
//       v < w - u
//
// Proof: By rearranging terms and substituting ">" for "<", we can assume
// that all values are non-negative.  Now clearly "w" is not the smallest
// value, so assume WLOG that "u" is the smallest.  We want to show that
// u < w - v in floating-point.  If v >= w/2, the calculation of w - v is
// exact since the result is smaller in magnitude than either input value,
// so the result holds.  Otherwise we have u <= v < w/2 and w - v >= w/2
// (even in floating point), so the result also holds.

// sum_equal reports whether u + v == w exactly.
fn sum_equal(u: f64, v: f64, w: f64) -> bool {
	(u+v == w) && (u == w-v) && (v == w-u)
}

#[derive(Clone,Copy,PartialEq,Debug)]
enum Axis {
    AxisU,
    AxisV,
}

impl Axis {
    fn value(self) -> u8 {
        match self {
            Axis::AxisU => 0,
            Axis::AxisV => 1,
        }
    }
}

type PointUVW = Point;

impl PointUVW {
    fn intersects_face(self) -> bool {
        // L intersects the [-1,1]x[-1,1] square in (u,v) if and only if the dot
        // products of N with the four corner vertices (-1,-1,1), (1,-1,1), (1,1,1),
        // and (-1,1,1) do not all have the same sign. This is true exactly when
        // |Nu| + |Nv| >= |Nw|. The code below evaluates this expression exactly.
        let u = self.0.x.abs();
        let v = self.0.y.abs();
        let w = self.0.z.abs();

        // We only need to consider the cases where u or v is the smallest value,
        // since if w is the smallest then both expressions below will have a
        // positive LHS and a negative RHS
        (v >= w-u) && (u >= w-v)
    }

    

    // intersects_opposite_edges reports whether a directed line L intersects two
    // opposite edges of a cube face F. This includs the case where L passes
    // exactly through a corner vertex of F. The directed line L is defined
    // by its normal N in the (u,v,w) coordinates of F.
    fn intersects_opposite_edges(self) -> bool {
        // The line L intersects opposite edges of the [-1,1]x[-1,1] (u,v) square if
        // and only exactly two of the corner vertices lie on each side of L. This
        // is true exactly when ||Nu| - |Nv|| >= |Nw|. The code below evaluates this
        // expression exactly.
        let u = self.0.x.abs();
        let v = self.0.y.abs();
        let w = self.0.z.abs();

        // If w is the smallest, the following line returns an exact result
        if (u-v).abs() != w {
            return (u-v).abs() >= w
        }

        // Otherwise u - v = w exactly, or w is not the smallest value. In either
        // case the following returns the correct result.
        if u >= v {
            return u-w >= v
        }
        v-w >= u
    }

    // exitAxis reports which axis the directed line L exits the cube face F on.
    // The directed line L is represented by its CCW normal N in the (u,v,w) coordinates
    // of F. It returns axisU if L exits through the u=-1 or u=+1 edge, and axisV if L exits
    // through the v=-1 or v=+1 edge. Either result is acceptable if L exits exactly
    // through a corner vertex of the cube face
    fn exit_axis(self) -> Axis {
        if self.intersects_opposite_edges() {
            // The line passes through through opposite edges of the face.
            // It exits through the v=+1 or v=-1 edge if the u-component of N has a
            // larger absolute magnitude than the v-component.
            if self.0.x.abs() >= self.0.y.abs(){
                return Axis::AxisV
            }
            return Axis::AxisU
        }

        // The line passes through through two adjacent edges of the face.
        // It exits the v=+1 or v=-1 edge if an even number of the components of N
        // are negative. We test this using signbit() rather than multiplication
        // to avoid the possibility of underflow.
        let signbit = |n: f64| -> u32 {
            if n.is_sign_negative() {
                1
            } else {
                0
            }
        };

        let x = signbit(self.0.x);
        let y = signbit(self.0.y);
        let z = signbit(self.0.z);
        if x^y^z == 0 {
            Axis::AxisV
        } else {
            Axis::AxisU
        }
    }

    // exitPoint returns the UV coordinates of the point where a directed line L (represented
    // by the CCW normal of this point), exits the cube face this point is derived from along
    // the given axis
    fn exit_point(self, a: Axis) -> r2::point::Point {
        if a == Axis::AxisU {
            let mut u = -1.0;
            if self.0.y > 0.0 {
                u = 1.0
            }
            return r2::point::Point{x: u, y: (-u*self.0.x - self.0.z) / self.0.y}
        }

        let mut v = -1.0;
        if self.0.x < 0.0 {
            v = 1.0
        }
        return r2::point::Point{x: (-v*self.0.y - self.0.z) / self.0.x, y: v}
    }
}

// clip_destination returns a score which is used to indicate if the clipped edge AB
// on the given face intersects the face at all. This function returns the score for
// the given endpoint, which is an integer ranging from 0 to 3. If the sum of the scores
// from both of the endpoints is 3 or more, then edge AB does not intersect this face.
//
// First, it clips the line segment AB to find the clipped destination B' on a given
// face. (The face is specified implicitly by expressing *all arguments* in the (u,v,w)
// coordinates of that face.) Second, it partially computes whether the segment AB
// intersects this face at all. The actual condition is fairly complicated, but it
// turns out that it can be expressed as a "score" that can be computed independently
// when clipping the two endpoints A and B.
fn clip_destination(a: PointUVW, b: PointUVW, scaled_n: PointUVW, a_tan: PointUVW, b_tan: PointUVW, scale_uv: f64) -> (r2::point::Point, i32) {
    // Optimization: if B is within the safe region of the face, use it.
    let max_save_uv_coord = 1.0 - FACE_CLIP_ERROR_UV_COORD;
    if b.0.z > 0.0 {
        let uv = r2::point::Point{x: b.0.x / b.0.z,  y: b.0.y / b.0.z };
        if uv.x.abs().max(uv.y.abs()) < max_save_uv_coord {
            return (uv, 0)
        } 
    }

    let mut uv = &(scaled_n.exit_point(scaled_n.exit_axis())) * scale_uv;
    let p: PointUVW = Point(r3::vector::Vector{x: uv.x, y: uv.y, z: 1.0});

    // Determine if the exit point B' is contained within the segment. We do this
	// by computing the dot products with two inward-facing tangent vectors at A
	// and B. If either dot product is negative, we say that B' is on the "wrong
	// side" of that point. As the point B' moves around the great circle AB past
	// the segment endpoint B, it is initially on the wrong side of B only; as it
	// moves further it is on the wrong side of both endpoints; and then it is on
	// the wrong side of A only. If the exit point B' is on the wrong side of
	// either endpoint, we can't use it; instead the segment is clipped at the
	// original endpoint B.
	//
	// We reject the segment if the sum of the scores of the two endpoints is 3
	// or more. Here is what that rule encodes:
	//  - If B' is on the wrong side of A, then the other clipped endpoint A'
	//    must be in the interior of AB (otherwise AB' would go the wrong way
	//    around the circle). There is a similar rule for A'.
	//  - If B' is on the wrong side of either endpoint (and therefore we must
	//    use the original endpoint B instead), then it must be possible to
	//    project B onto this face (i.e., its w-coordinate must be positive).
	//    This rule is only necessary to handle certain zero-length edges (A=B).
    let mut score = 0;
    if (p.0 - a.0).dot(&a_tan.0) < 0.0 {
        score = 2; // B' is on wrong side of A.
    } else if (p.0 - b.0).dot(&b_tan.0) < 0.0{
        score = 1; // B' is on wrong side of B.
    }

    if score > 0 { // B' is not in the interior of AB.
        if b.0.z <= 0.0 {
            score = 3; // B cannot be projected onto this face.
        } else {
            let v = b.0;
            uv = r2::point::Point{x: v.x / v.z, y: v.y / v.z}
        }
    }

    (uv, score)
}

// update_endpoint returns the interval with the specified endpoint updated to
// the given value. If the value lies beyond the opposite endpoint, nothing is
// changed and false is returned
fn update_endpoint(mut bound: r1::interval::Interval, high_endpoint: bool, value: f64) -> (r1::interval::Interval, bool) {
    if !high_endpoint {
        if bound.hi < value {
            return (bound, false)
        }
        if bound.lo < value {
            bound.lo = value
        }
        return (bound, true)
    }

    if bound.lo > value {
        return (bound, false)
    }
    if bound.hi > value {
        bound.hi = value
    }
    (bound, true)
}

// clip_bound_axis returns the clipped versions of the bounding intervals for the given
// axes for the line segment from (a0,a1) to (b0,b1) so that neither extends beyond the
// given clip interval. negSlope is a precomputed helper variable that indicates which
// diagonal of the bounding box is spanned by AB; it is false if AB has positive slope,
// and true if AB has negative slope. If the clipping interval doesn't overlap the bounds,
// false is returned
fn clip_bound_axis(a0: f64, b0: f64, mut bound0: r1::interval::Interval, a1: f64, b1: f64, bound1: r1::interval::Interval, neg_slope: bool, clip: r1::interval::Interval) 
    -> (r1::interval::Interval, r1::interval::Interval, bool) {

    if bound0.lo < clip.lo {
        // If the upper bound is below the clips lower bound, there is nothing to do.
        if bound0.hi < clip.lo {
            return (bound0, bound1, false)
        }
        // narrow the intervals lower bound to the clip bound
        bound0.lo = clip.lo;
        let (bound1, updated) = update_endpoint(bound1, neg_slope, interpolate(clip.lo, a0, b0, a1, b1));
        if !updated {
            return (bound0, bound1, false) 
        }
    }

    if bound0.hi > clip.hi {
        // If the lower bound is above the clips upper bound, there is nothing to do.
        if bound0.lo > clip.hi {
            return (bound0, bound1, false)
        }
        // narrow the intervals upper bound to the clip bound.
        bound0.hi = clip.hi;
        let (bound1, updated) = update_endpoint(bound1, !neg_slope, interpolate(clip.hi, a0, b0, a1, b1));
        if !updated {
            return (bound0, bound1, false)
        }
    }
    return (bound0, bound1, true)
}

// edge_intersects_rect reports whether the edge defined by AB intersects the
// given closed rectangle to within the error bound.
#[allow(dead_code)]
fn edge_intersects_rec(a: r2::point::Point, b: r2::point::Point, r: &r2::rect::Rect) -> bool {
    // First check whether the bounds of a Rect around AB intersects the given rect.
    if !r.intersects(&(r2::rect::Rect::from_points(&[a,b]))) {
        return false
    }

    // Otherwise AB intersects the rect if and only if all four vertices of rect
	// do not lie on the same side of the extended line AB. We test this by finding
	// the two vertices of rect with minimum and maximum projections onto the normal
	// of AB, and computing their dot products with the edge normal.
    let n = (b - a).ortho();
    let (mut i, mut j):(isize, isize) = (0,0);
    if n.x >= 0.0 {
        i = 1;
    }
    if n.y >= 0.0 {
        j = 1;
    }
    let max = n.dot(&(r.vertex_ij(i,j) - a));
    let min = n.dot(&(r.vertex_ij(1-i,1-j) - a));

    return (max >= 0.0) && (min <= 0.0)
}

// clippedEdgeBound returns the bounding rectangle of the portion of the edge defined
// by AB intersected by clip. The resulting bound may be empty. This is a convenience
// function built on top of clipEdgeBound.
#[allow(dead_code)]
fn clipped_edge_bound(a: r2::point::Point, b: r2::point::Point, clip: r2::rect::Rect) -> r2::rect::Rect  {
    let bound = r2::rect::Rect::from_points(&[a,b]);
    let (b1, intersects) = clip_edge_bound(a, b, clip, bound);
    if intersects {
        b1
    } else {
        r2::rect::EMPTY
    }
}

// clip_edge_bound clips an edge AB to sequence of rectangles efficiently.
// It represents the clipped edges by their bounding boxes rather than as a pair of
// endpoints. Specifically, let A'B' be some portion of an edge AB, and let bound be
// a tight bound of A'B'. This function returns the bound that is a tight bound
// of A'B' intersected with a given rectangle. If A'B' does not intersect clip,
// it returns false and the original bound.
#[allow(dead_code)]
fn clip_edge_bound(a: r2::point::Point, b: r2::point::Point, clip: r2::rect::Rect, bound: r2::rect::Rect) -> (r2::rect::Rect, bool) {
	// negSlope indicates which diagonal of the bounding box is spanned by AB: it
	// is false if AB has positive slope, and true if AB has negative slope. This is
	// used to determine which interval endpoints need to be updated each time
	// the edge is clipped
    let neg_slope = (a.x > b.x) != (a.y > b.y);
    let (b0_x, b0_y, up1) = clip_bound_axis(a.x, b.x, bound.x, a.y, b.y, bound.y, neg_slope, clip.x);
    if !up1 {
        return (bound, false)
    }
    let (b1_x, b1_y, up2) = clip_bound_axis(a.x, b.x, bound.x, a.y, b.y, bound.y, neg_slope, clip.x);
    if !up2 {
        return (r2::rect::Rect{x: b0_x, y: b0_y}, false)
    }
    return (r2::rect::Rect{x: b1_x, y: b1_y}, true)
}

// interpolate returns a value with the same combination of a1 and b1 as the
// given value x is of a and b. This function makes the following guarantees:
//  - If x == a, then x1 = a1 (exactly).
//  - If x == b, then x1 = b1 (exactly).
//  - If a <= x <= b, then a1 <= x1 <= b1 (even if a1 == b1).
// This requires a != b.
fn interpolate(x: f64, a: f64, b: f64, a1: f64, b1: f64) -> f64 {
    if (a-x).abs() <= (b-x).abs() {
        return a1 + (b1-a1) * (x-a) / (b-a)
    }
    b1 + (a1-b1) * (x-b) / (a-b)
}

// FaceSegment represents an edge AB clipped to an S2 cube face. It is
// represented by a face index and a pair of (u,v) coordinates.
#[derive(Clone,Copy,PartialEq,Debug)]
pub struct FaceSegment {
    pub face: u8,
    pub a: r2::point::Point,
    pub b: r2::point::Point,
}

// face_segments subdivides the given edge AB at every point where it crosses the
// boundary between two S2 cube faces and returns the corresponding FaceSegments.
// The segments are returned in order from A toward B. The input points must be
// unit length.
//
// This function guarantees that the returned segments form a continuous path
// from A to B, and that all vertices are within faceClipErrorUVDist of the
// line AB. All vertices lie within the [-1,1]x[-1,1] cube face rectangles.
// The results are consistent with Sign, i.e. the edge is well-defined even its
// endpoints are antipodal.
// TODO(roberts): Extend the implementation of PointCross so that this is true.
pub fn face_segments(a: Point, b: Point) -> Vec<FaceSegment> {
    // Fast path: both endpoints are on the same face
    let (a_face, a_x, a_y) = stuv::xyz_to_face_uv(&a.0);
    let (b_face, b_x, b_y) = stuv::xyz_to_face_uv(&b.0);

    if a_face == b_face {
        return vec![
            FaceSegment {
                face: a_face,
                a: r2::point::Point{ x: a_x, y: a_y },
                b: r2::point::Point{ x: b_x, y: b_y },
            }
        ]
    }

	// Starting at A, we follow AB from face to face until we reach the face
	// containing B. The following code is designed to ensure that we always
	// reach B, even in the presence of numerical errors.
	//
	// First we compute the normal to the plane containing A and B. This normal
	// becomes the ultimate definition of the line AB; it is used to resolve all
	// questions regarding where exactly the line goes. Unfortunately due to
	// numerical errors, the line may not quite intersect the faces containing
	// the original endpoints. We handle this by moving A and/or B slightly if
	// necessary so that they are on faces intersected by the line AB.
    let mut segment = FaceSegment {
        face: a_face,
        a: r2::point::Point{ x: a_x, y: a_y },
        b: r2::point::Point{ x: b_x, y: b_y },
    };

    #[allow(unused)]
    let ab = a.cross(&b);

    let (a_face, segment_a) = move_origin_to_valid_face(a_face, a, ab, segment.a);
    let (b_face, segment_b) = move_origin_to_valid_face(b_face, b, ab * -1.0, segment.b);

    segment = FaceSegment{
        face: a_face,
        a: segment_a,
        b: segment_b,
    };

    let mut segments = Vec::new();
    let b_saved = segment_b;
    let mut face = a_face;
    while face != b_face {
        // Complete the current segment by finding the point where AB
		// exits the current face

        let n: PointUVW = stuv::face_xyz_to_uvw(face, &ab);
        let exit_axis = n.exit_axis();
        segment.b = n.exit_point(exit_axis);
        segments.push(segment.clone());
        
        // Compute the next face intersected by AB, and translate the exit
		// point of the current segment into the (u,v) coordinates of the
		// next face. This becomes the first point of the next segment. 
        let exit_xyz = stuv::face_uv_to_xyz(face, segment.b.x, segment.b.y);
        face = next_face(face, segment.b, exit_axis, n, b_face);
        let exit_uvw = stuv::face_xyz_to_uvw(face, &Point(exit_xyz));
        segment.face = face;
        segment.a = r2::point::Point{x: exit_uvw.0.x, y: exit_uvw.0.y};
    }
    segment.b = b_saved;
    segments.push(segment);
    segments
}


// moveOriginToValidFace updates the origin point to a valid face if necessary.
// Given a line segment AB whose origin A has been projected onto a given cube
// face, determine whether it is necessary to project A onto a different face
// instead. This can happen because the normal of the line AB is not computed
// exactly, so that the line AB (defined as the set of points perpendicular to
// the normal) may not intersect the cube face containing A. Even if it does
// intersect the face, the exit point of the line from that face may be on
// the wrong side of A (i.e., in the direction away from B). If this happens,
// we reproject A onto the adjacent face where the line AB approaches A most
// closely. This moves the origin by a small amount, but never more than the
// error tolerances.
fn move_origin_to_valid_face(mut face: u8, a: Point, ab: Point, mut a_uv: r2::point::Point) -> (u8, r2::point::Point) {
    // Fast path: if the origin is sufficiently far inside the face, it is
	// always safe to use it.
    let max_save_uv_coord = 1.0 - FACE_CLIP_ERROR_UV_COORD;
    if a_uv.x.abs().max(a_uv.y.abs()) <= max_save_uv_coord {
        return (face, a_uv)
    }
    
    // Otherwise check whether the normal AB even intersects this face
    let n: PointUVW = stuv::face_xyz_to_uvw(face, &ab);
    if n.intersects_face() {
		// Check whether the point where the line AB exits this face is on the
		// wrong side of A (by more than the acceptable error tolerance).
        let uv = n.exit_point(n.exit_axis());
        let exit = stuv::face_uv_to_xyz(face, uv.x, uv.y);
        let a_tan = ab.normalize().cross(&a);

		// We can use the given face
        if (exit - a.0).dot(&a_tan.0) >= -FACE_CLIP_ERROR_RADIANS {
            return (face, a_uv)
        }
    }

    // Otherwise we reproject A to the nearest adjacent face. (If line AB does
	// not pass through a given face, it must pass through all adjacent faces.)
    let mut dir: u8 = 0;
    if a_uv.x.abs() >= a_uv.y.abs(){
        // U-axis
        if a_uv.x > 0.0 {
            dir =1;
        }
        face = stuv::uvw_face(face, 0, dir);
    } else {
        // V-axis
        if a_uv.y > 0.0 {
            dir = 1;
        }
        face = stuv::uvw_face(face, 1, dir);
    }

    let (a_uv_x, a_uv_y) = stuv::valid_face_xyz_to_uv(face, &a.0);
    a_uv.x = a_uv_x;
    a_uv.y = a_uv_y;
    a_uv.x = -1.0f64.max(1.0f64.min(a_uv.x));
    a_uv.y = -1.0f64.max(1.0f64.min(a_uv.y));
    (face, a_uv)
}

// nextFace returns the next face that should be visited by FaceSegments, given that
// we have just visited face and we are following the line AB (represented
// by its normal N in the (u,v,w) coordinates of that face). The other
// arguments include the point where AB exits face, the corresponding
// exit axis, and the target face containing the destination point B
fn next_face(face: u8, exit: r2::point::Point, axis: Axis, n: PointUVW, target_face: u8) -> u8 {
    // this bit is to work around C++ cleverly casting bools to ints for you.
    let mut exit_a = exit.x;
    let mut exit_1minus_a = exit.y;

    if axis == Axis::AxisV {
        exit_a = exit.y;
        exit_1minus_a = exit.x;
    }

    let exit_a_pos = if exit_a > 0.0 {
        1
    } else {
        0
    };

    let exit_1minus_a_pos = if exit_1minus_a > 0.0 {
        1
    } else {
        0
    };

    // We return the face that is adjacent to the exit point along the given
	// axis. If line AB exits *exactly* through a corner of the face, there are
	// two possible next faces. If one is the target face containing B, then
	// we guarantee that we advance to that face directly.
	//
	// The three conditions below check that (1) AB exits approximately through
	// a corner, (2) the adjacent face along the non-exit axis is the target
	// face, and (3) AB exits *exactly* through the corner. (The sumEqual
	// code checks whether the dot product of (u,v,1) and n is exactly zero.)
    if exit_1minus_a.abs() == 1.0 && stuv::uvw_face(face, 1-axis.value(), exit_1minus_a_pos) == target_face && sum_equal(exit.x*n.0.x, exit.y*n.0.y, -n.0.z) {
        return target_face
    }

    // Otherwise return the face that is adjacent to the exit point in the
	// direction of the exit axis
    stuv::uvw_face(face, axis.value(), exit_a_pos)
}


#[cfg(test)]
pub mod test {
    use super::*;
    use point::Point;
    use r3;

    #[test]
    fn test_edge_clipping_intersects_face() {
        /*
		{pointUVW{r3.Vector{-3.88578e-16, -math.Sqrt(2.0 / 3.0), math.Sqrt(2.0 / 3.0)}}, true}
        */
        assert!((Point(r3::vector::Vector{x: 2.05335e-06, y: 3.91604e-22, z: 2.90553e-06}) as PointUVW).intersects_face() == false);
        assert!((Point(r3::vector::Vector{x: -3.91604e-22, y: -2.05335e-06, z: -2.90553e-06}) as PointUVW).intersects_face() == false);
        assert!((Point(r3::vector::Vector{x: 0.169258, y: -0.169258, z: 0.664013}) as PointUVW).intersects_face() == false);
        assert!((Point(r3::vector::Vector{x: 0.169258, y: -0.169258, z: 0.664013}) as PointUVW).intersects_face() == false);
        assert!((Point(r3::vector::Vector{x: (2.0_f64 / 3.0_f64).sqrt() , y: -(2.0_f64 / 3.0_f64).sqrt(), z: 3.88578e-16}) as PointUVW).intersects_face() == true);
        assert!((Point(r3::vector::Vector{x: 3.88578e-16 , y: -(2.0_f64 / 3.0_f64).sqrt(), z: (2.0_f64 / 3.0_f64).sqrt()}) as PointUVW).intersects_face() == true);
    }
    
    #[test]
    fn test_edge_clipping_intersects_opposite_edges() {
        assert!((Point(r3::vector::Vector{x: 0.169258, y: -0.169258, z: 0.664013}) as PointUVW).intersects_opposite_edges() == false);
        assert!((Point(r3::vector::Vector{x: 0.169258, y: -0.169258, z: -0.664013}) as PointUVW).intersects_opposite_edges() == false);
        assert!((Point(r3::vector::Vector{x: (4.0_f64 / 3.0_f64).sqrt(), y: 0_f64, z: -(4_f64 / 3_f64).sqrt()}) as PointUVW).intersects_opposite_edges() == true);
        assert!((Point(r3::vector::Vector{x: (4.0_f64 / 3.0_f64).sqrt(), y: 0_f64, z: (4_f64 / 3_f64).sqrt()}) as PointUVW).intersects_opposite_edges() == true);
        assert!((Point(r3::vector::Vector{x: -(2.0_f64 / 3.0_f64).sqrt(), y: -(2.0_f64 / 3.0_f64).sqrt(), z: 1.66533453694e-16}) as PointUVW).intersects_opposite_edges() == false);
        assert!((Point(r3::vector::Vector{x: (2.0_f64 / 3.0_f64).sqrt(), y: (2.0_f64 / 3.0_f64).sqrt(), z: -1.66533453694e-16}) as PointUVW).intersects_opposite_edges() == false);
    }

    #[test]
    fn test_edge_clipping_exit_point() {

    }
}
