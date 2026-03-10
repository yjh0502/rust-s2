// Copyright 2017 Google Inc. All rights reserved.
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

use std::cmp::Ordering;
use std::ops::{Add, Mul, Sub};

use crate::consts::{DBL_EPSILON, EPSILON};
use crate::point::{ordered_ccw, Point};
use crate::r3::precisevector::PreciseVector;
use crate::r3::vector::*;
pub(crate) use crate::s2::edge_crosser::EdgeCrosser;

// intersectionError can be set somewhat arbitrarily, because the algorithm
// uses more precision if necessary in order to achieve the specified error.
// The only strict requirement is that intersectionError >= dblEpsilon
// radians. However, using a larger error tolerance makes the algorithm more
// efficient because it reduces the number of cases where exact arithmetic is
// needed.
static INTERSECTION_ERROR: f64 = 8.0 * DBL_EPSILON;

// intersectionMergeRadius is used to ensure that intersection points that
// are supposed to be coincident are merged back together into a single
// vertex. This is required in order for various polygon operations (union,
// intersection, etc) to work correctly. It is twice the intersection error
// because two coincident intersection points might have errors in
// opposite directions.
static INTERSECTION_MERGE_RADIUS: f64 = 16.0 * DBL_EPSILON;

// A Crossing indicates how edges cross.
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Crossing {
    Cross,
    Maybe,
    DoNotCross,
}

impl std::fmt::Display for Crossing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Crossing::Cross => write!(f, "Cross"),
            Crossing::Maybe => write!(f, "MaybeCross"),
            Crossing::DoNotCross => write!(f, "DoNotCross"),
        }
    }
}

impl Crossing {
    // CrossingSign reports whether the edge AB intersects the edge CD.
    // If AB crosses CD at a point that is interior to both edges, Cross is returned.
    // If any two vertices from different edges are the same it returns MaybeCross.
    // Otherwise it returns DoNotCross.
    // If either edge is degenerate (A == B or C == D), the return value is MaybeCross
    // if two vertices from different edges are the same and DoNotCross otherwise.
    //
    // Properties of CrossingSign:
    //
    //	(1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
    //	(2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
    //	(3) CrossingSign(a,b,c,d) == MaybeCross if a==c, a==d, b==c, b==d
    //	(3) CrossingSign(a,b,c,d) == DoNotCross or MaybeCross if a==b or c==d
    //
    // This method implements an exact, consistent perturbation model such
    // that no three points are ever considered to be collinear. This means
    // that even if you have 4 points A, B, C, D that lie exactly in a line
    // (say, around the equator), C and D will be treated as being slightly to
    // one side or the other of AB. This is done in a way such that the
    // results are always consistent (see RobustSign).
    pub fn crossing_sign(a: &Point, b: &Point, c: &Point, d: &Point) -> Crossing {
        let mut crosser = EdgeCrosser::new_chain_edge_crosser(a, b, c);
        crosser.chain_crossing_sign(d)
    }
}

// VertexCrossing reports whether two edges "cross" in such a way that point-in-polygon
// containment tests can be implemented by counting the number of edge crossings.
//
// Given two edges AB and CD where at least two vertices are identical
// (i.e. CrossingSign(a,b,c,d) == 0), the basic rule is that a "crossing"
// occurs if AB is encountered after CD during a CCW sweep around the shared
// vertex starting from a fixed reference point.
//
// Note that according to this rule, if AB crosses CD then in general CD
// does not cross AB. However, this leads to the correct result when
// counting polygon edge crossings. For example, suppose that A,B,C are
// three consecutive vertices of a CCW polygon. If we now consider the edge
// crossings of a segment BP as P sweeps around B, the crossing number
// changes parity exactly when BP crosses BA or BC.
//
// Useful properties of VertexCrossing (VC):
//
//	(1) VC(a,a,c,d) == VC(a,b,c,c) == false
//	(2) VC(a,b,a,b) == VC(a,b,b,a) == true
//	(3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
//	(3) If exactly one of a,b equals one of c,d, then exactly one of
//	    VC(a,b,c,d) and VC(c,d,a,b) is true
//
// It is an error to call this method with 4 distinct vertices.
pub fn vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> bool {
    // If A == B or C == D there is no intersection. We need to check this
    // case first in case 3 or more input points are identical.
    if a == b || c == d {
        return false;
    }

    // If any other pair of vertices is equal, there is a crossing if and only
    // if ordered_ccw indicates that the edge AB is further CCW around the
    // shared vertex O (either A or B) than the edge CD, starting from an
    // arbitrary fixed reference point.

    // Optimization: if AB=CD or AB=DC, we can avoid most of the calculations.
    if a == c {
        return (b == d) || ordered_ccw(&a.reference_dir(), d, b, a);
    } else if b == d {
        return ordered_ccw(&b.reference_dir(), c, a, b);
    } else if a == d {
        return (b == c) || ordered_ccw(&a.reference_dir(), c, b, a);
    } else if b == c {
        return ordered_ccw(&b.reference_dir(), d, a, b);
    }

    false
}

// EdgeOrVertexCrossing is a convenience function that calls CrossingSign to
// handle cases where all four vertices are distinct, and VertexCrossing to
// handle cases where two or more vertices are the same. This defines a crossing
// function such that point-in-polygon containment tests can be implemented
// by simply counting edge crossings.
pub fn edge_or_vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> bool {
    match Crossing::crossing_sign(a, b, c, d) {
        Crossing::DoNotCross => false,
        Crossing::Cross => true,
        Crossing::Maybe => vertex_crossing(a, b, c, d),
    }
}

// Intersection returns the intersection point of two edges AB and CD that cross
// (CrossingSign(a,b,c,d) == Crossing).
//
// Useful properties of Intersection:
//
//	(1) Intersection(b,a,c,d) == Intersection(a,b,d,c) == Intersection(a,b,c,d)
//	(2) Intersection(c,d,a,b) == Intersection(a,b,c,d)
//
// The returned intersection point X is guaranteed to be very close to the
// true intersection point of AB and CD, even if the edges intersect at a
// very small angle.
pub fn intersection(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Point {
    // It is difficult to compute the intersection point of two edges accurately
    // when the angle between the edges is very small. Previously we handled
    // this by only guaranteeing that the returned intersection point is within
    // intersectionError of each edge. However, this means that when the edges
    // cross at a very small angle, the computed result may be very far from the
    // true intersection point.
    //
    // Instead this fntion now guarantees that the result is always within
    // intersectionError of the true intersection. This requires using more
    // sophisticated techniques and in some cases extended precision.
    //
    //  - intersectionStable computes the intersection point using
    //    projection and interpolation, taking care to minimize cancellation
    //    error.
    //
    //  - intersectionExact computes the intersection point using precision
    //    arithmetic and converts the final result back to an Point.
    let mut pt =
        intersection_stable(a0, a1, b0, b1).unwrap_or_else(|_e| intersection_exact(a0, a1, b0, b1));

    // Make sure the intersection point is on the correct side of the sphere.
    // Since all vertices are unit length, and edges are less than 180 degrees,
    // (a0 + a1) and (b0 + b1) both have positive dot product with the
    // intersection point.  We use the sum of all vertices to make sure that the
    // result is unchanged when the edges are swapped or reversed.
    if pt.0.dot(&a0.0.add(a1.0.add(b0.0.add(b1.0)))) < 0.0 {
        pt = Point(pt.0.mul(-1));
    }

    pt
}

// Computes the cross product of two vectors, normalized to be unit length.
// Also returns the length of the cross
// product before normalization, which is useful for estimating the amount of
// error in the result.  For numerical stability, the vectors should both be
// approximately unit length.
fn robust_normal_with_length(x: &Vector, y: &Vector) -> (Vector, f64) {
    // This computes 2 * (x.cross(y)), but has much better numerical
    // stability when x and y are unit length.
    let tmp = x.sub(y).cross(&x.add(y));
    let length = tmp.norm();
    if length != 0f64 {
        (tmp.mul(1.0 / length), 0.5 * length)
    } else {
        (Vector::default(), 0.5 * length)
    }
}

/*
// intersectionSimple is not used by the C++ so it is skipped here.
*/

// projection returns the projection of a_norm onto X (x.dot(a_norm)), and a bound
// on the error in the result. a_norm is not necessarily unit length.
//
// The remaining parameters (the length of a_norm (a_norm_len) and the edge endpoints
// a0 and a1) allow this dot product to be computed more accurately and efficiently.
fn projection(x: &Vector, a_norm: &Vector, a_norm_len: f64, a0: &Point, a1: &Point) -> (f64, f64) {
    // The error in the dot product is proportional to the lengths of the input
    // vectors, so rather than using x itself (a unit-length vector) we use
    // the vectors from x to the closer of the two edge endpoints. This
    // typically reduces the error by a huge factor.
    let x0 = x.sub(&a0.0);
    let x1 = x.sub(&a1.0);
    let x0_dist2 = x0.norm2();
    let x1_dist2 = x1.norm2();

    // If both distances are the same, we need to be careful to choose one
    // endpoint deterministically so that the result does not change if the
    // order of the endpoints is reversed.
    let mut dist = 0.0;
    let mut proj = 0.0;

    if x0_dist2 < x1_dist2 || (x0_dist2 == x1_dist2 && x0.cmp(&x1) == Ordering::Less) {
        dist = x0_dist2.sqrt();
        proj = x0.dot(&a_norm);
    } else {
        dist = x1_dist2.sqrt();
        proj = x1.dot(&a_norm)
    }

    // This calculation bounds the error from all sources: the computation of
    // the normal, the subtraction of one endpoint, and the dot product itself.
    // dblError appears because the input points are assumed to be
    // normalized in double precision.
    //
    // For reference, the bounds that went into this calculation are:
    // ||N'-N|| <= ((1 + 2 * sqrt(3))||N|| + 32 * sqrt(3) * dblError) * epsilon
    // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * epsilon
    // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * epsilon
    let bound = (((3.5 + 2.0 * 3.0_f64.sqrt()) * a_norm_len + 32. * 3.0_f64.sqrt() * DBL_EPSILON)
        * dist
        + 1.5 * proj.abs())
        * EPSILON;

    return (proj, bound);
}

// compareEdges reports whether (a0,a1) is less than (b0,b1) with respect to a total
// ordering on edges that is invariant under edge reversals.
fn compare_edges(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> bool {
    let (a0, _a1) = match a0.0.cmp(&a1.0) {
        Ordering::Less => (a0, a1),
        Ordering::Equal => (a1, a0),
        Ordering::Greater => (a1, a0),
    };

    let (b0, b1) = match b0.0.cmp(&b1.0) {
        Ordering::Less => (b0, b1),
        Ordering::Equal => (b1, b0),
        Ordering::Greater => (b1, b0),
    };

    a0.0 < b0.0 || (a0 == b0 && b0.0 < b1.0)
}

// intersectionStable returns the intersection point of the edges (a0,a1) and
// (b0,b1) if it can be computed to within an error of at most intersectionError
// by this fntion.
//
// The intersection point is not guaranteed to have the correct sign because we
// choose to use the longest of the two edges first. The sign is corrected by
// Intersection.
fn intersection_stable(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Result<Point, Point> {
    // Sort the two edges so that (a0,a1) is longer, breaking ties in a
    // deterministic way that does not depend on the ordering of the endpoints.
    // This is desirable for two reasons:
    //  - So that the result doesn't change when edges are swapped or reversed.
    //  - It reduces error, since the first edge is used to compute the edge
    //    normal (where a longer edge means less error), and the second edge
    //    is used for interpolation (where a shorter edge means less error).
    let a_len2 = a1.0.sub(a0.0).norm2();
    let b_len2 = b1.0.sub(b0.0).norm2();
    if a_len2 < b_len2 || (a_len2 == b_len2 && compare_edges(a0, a1, b0, b1)) {
        return intersect_stable_sorted(b0, b1, a0, a1);
    }

    intersect_stable_sorted(a0, a1, b0, b1)
}

// intersectionStableSorted is a helper fntion for intersectionStable.
// It expects that the edges (a0,a1) and (b0,b1) have been sorted so that
// the first edge passed in is longer.
fn intersect_stable_sorted(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Result<Point, Point> {
    let pt = Point::default();

    // Compute the normal of the plane through (a0, a1) in a stable way.
    let a_norm = a0.0.sub(a1.0).cross(&a0.0.add(a1.0));
    let a_norm_len = a_norm.norm();
    let b_len = b1.sub(b0).norm();

    // Compute the projection (i.e., signed distance) of b0 and b1 onto the
    // plane through (a0, a1).  Distances are scaled by the length of a_norm.
    let (b0_dist, b0_error) = projection(&b0.0, &a_norm, a_norm_len, &a0, &a1);
    let (b1_dist, b1_error) = projection(&b1.0, &a_norm, a_norm_len, &a0, &a1);

    // The total distance from b0 to b1 measured perpendicularly to (a0,a1) is
    // |b0Dist - b1Dist|.  Note that b0Dist and b1Dist generally have
    // opposite signs because b0 and b1 are on opposite sides of (a0, a1).  The
    // code below finds the intersection point by interpolating along the edge
    // (b0, b1) to a fractional distance of b0Dist / (b0Dist - b1Dist).
    //
    // It can be shown that the maximum error in the interpolation fraction is
    //
    //   (b0Dist * b1Error - b1Dist * b0Error) / (dist_sum * (dist_sum - error_sum))
    //
    // We save ourselves some work by scaling the result and the error bound by
    // "dist_sum", since the result is normalized to be unit length anyway.
    let dist_sum = (b0_dist - b1_dist).abs();
    let error_sum = b0_error + b1_error;
    if dist_sum <= error_sum {
        return Err(pt); // Error is unbounded in this case.
    }

    let x = b1.mul(b0_dist).sub(b0.mul(b1_dist));
    let err = b_len * (b0_dist * b1_error - b1_dist * b0_error).abs() / (dist_sum - error_sum)
        + 2.0 * dist_sum * DBL_EPSILON;

    // Finally we normalize the result, compute the corresponding error, and
    // check whether the total error is acceptable.
    let x_len = x.norm();
    if err > (INTERSECTION_ERROR - EPSILON) * x_len {
        return Err(pt);
    }

    Ok(x.mul(1.0 / x_len))
}

// intersectionExact returns the intersection point of (a0, a1) and (b0, b1)
// using precise arithmetic. Note that the result is not exact because it is
// rounded down to double precision at the end. Also, the intersection point
// is not guaranteed to have the correct sign (i.e., the return value may need
// to be negated).
fn intersection_exact(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Point {
    // Since we are using presice arithmetic, we don't need to worry about
    // numerical stability.
    let a0_p: PreciseVector = a0.0.into();
    let a1_p: PreciseVector = a1.0.into();
    let b0_p: PreciseVector = b0.0.into();
    let b1_p: PreciseVector = b1.0.into();
    let a_norm_p = a0_p.cross(&a1_p);
    let b_norm_p = b0_p.cross(&b1_p);
    let x_p = a_norm_p.cross(&b_norm_p);

    // The final Normalize() call is done in double precision, which creates a
    // directional error of up to 2*dblError. (Precise conversion and Normalize()
    // each contribute up to dblError of directional error.)
    let mut x: Vector = x_p.into();

    if x == Vector::default() {
        // The two edges are exactly collinear, but we still consider them to be
        // "crossing" because of simulation of simplicity. Out of the four
        // endpoints, exactly two lie in the interior of the other edge. Of
        // those two we return the one that is lexicographically smallest.
        x = Vector::new(10.0, 10.0, 10.0); // Greater than any valid S2Point

        let a_norm = Point(a_norm_p.into());
        let b_norm = Point(b_norm_p.into());

        if ordered_ccw(b0, a0, b1, &b_norm) && a0.0 < x {
            return *a0;
        }
        if ordered_ccw(b0, a1, b1, &b_norm) && a1.0 < x {
            return *a1;
        }
        if ordered_ccw(a0, b0, a1, &a_norm) && b0.0 < x {
            return *b0;
        }
        if ordered_ccw(a0, b1, a1, &a_norm) && b1.0 < x {
            return *b1;
        }
    }

    Point(x)
}

// AngleContainsVertex reports if the angle ABC contains its vertex B.
// Containment is defined such that if several polygons tile the region around
// a vertex, then exactly one of those polygons contains that vertex.
// Returns false for degenerate angles of the form ABA.
//
// Note that this method is not sufficient to determine vertex containment in
// polygons with duplicate vertices (such as the polygon ABCADE).  Use
// ContainsVertexQuery for such polygons. AngleContainsVertex(a, b, c)
// is equivalent to using ContainsVertexQuery as follows:
//
//	ContainsVertexQuery query(b);
//	query.addEdge(a, -1);  // incoming
//	query.addEdge(c, 1);   // outgoing
//	return query.ContainsVertex() > 0;
//
// Useful properties of AngleContainsVertex:
//
//	(1) AngleContainsVertex(a,b,a) == false
//	(2) AngleContainsVertex(a,b,c) == !AngleContainsVertex(c,b,a) unless a == c
//	(3) Given vertices v_1 ... v_k ordered cyclically CCW around vertex b,
//	    AngleContainsVertex(v_{i+1}, b, v_i) is true for exactly one value of i.
//
// REQUIRES: a != b && b != c
pub fn angle_contains_vertex(a: &Point, b: &Point, c: &Point) -> bool {
    // A loop with consecutive vertices A, B, C contains vertex B if and only if
    // the fixed vector R = referenceDir(B) is contained by the wedge ABC.  The
    // wedge is closed at A and open at C, i.e. the point B is inside the loop
    // if A = R but not if C = R.
    //
    // Note that the test below is written so as to get correct results when the
    // angle ABC is degenerate. If A = C or C = R it returns false, and
    // otherwise if A = R it returns true.
    return !ordered_ccw(&b.reference_dir(), c, a, b);
}

// TODO(roberts): Differences from C++
// fn RobustCrossProd(a, b Point) Point
// fn symbolicCrossProd(a, b Point) Point
// fn exactCrossProd(a, b Point) Point
// fn SignedVertexCrossing(a, b, c, d Point) int
// fn isNormalizable(p Point) bool
// fn ensureNormalizable(p Point) Point
// fn normalizableFromPrecise(p r3.PreciseVector) Point

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::DBL_EPSILON;
    use crate::r3::vector::Vector;
    use crate::s1::angle::Angle;

    // Helper function to create a Point from coordinates
    fn point_from_coords(x: f64, y: f64, z: f64) -> Point {
        Point(Vector::new(x, y, z).normalize())
    }

    // Helper function for the tests to return a positively
    // oriented intersection Point of the two line segments (a0,a1) and (b0,b1).
    fn test_intersection_exact(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Point {
        let mut x = intersection_exact(a0, a1, b0, b1);
        if x.0.dot(&a0.0.add(a1.0.add(b0.0.add(b1.0)))) < 0.0 {
            x = Point(x.0.mul(-1.0));
        }
        x
    }

    // DistanceFromSegment returns the distance from point X to line segment AB.
    fn distance_from_segment(x: &Point, a: &Point, b: &Point) -> Angle {
        // We compute the distance using the standard formula for the distance from a
        // point X to a line segment AB:
        //
        // (1) If the line AB is degenerate, return the distance to the closest endpoint.
        //
        // (2) Otherwise, if the projection of X onto the line AB is outside the range
        // [0,1], then return the distance to the closest endpoint.
        //
        // (3) Otherwise, return the perpendicular distance from X to the line AB.
        //
        // Note that the projection parameter is the dot product of the unit vector in
        // the direction AB with the vector AX, where both vectors are treated as
        // Euclidean (not spherical).

        println!(
            "Computing distance from {:?} to segment [{:?}, {:?}]",
            x, a, b
        );

        if a == b {
            let dist = x.distance(a);
            println!(
                "Degenerate segment, returning distance to endpoint: {}",
                dist.0
            );
            return dist;
        }

        let ab = b.0.sub(a.0);
        let ax = x.0.sub(a.0);
        let ab_norm2 = ab.norm2();
        let ax_dot_ab = ax.dot(&ab);

        println!("ab_norm2: {}", ab_norm2);
        println!("ax_dot_ab: {}", ax_dot_ab);

        // Handle cases (2) and (3).
        if ax_dot_ab <= 0.0 {
            let dist = x.distance(a);
            println!("Closest to endpoint A, distance: {}", dist.0);
            return dist;
        }
        if ax_dot_ab >= ab_norm2 {
            let dist = x.distance(b);
            println!("Closest to endpoint B, distance: {}", dist.0);
            return dist;
        }

        // The closest point is the projection of X onto AB.
        // For numerical stability, we need to be careful with the normalization step
        // when dealing with points that are extremely close to each other.
        let p = a.0.add(ab.mul(ax_dot_ab / ab_norm2));

        // Check if p is already very close to a unit vector
        let p_norm = p.norm();
        let normalized_p = if (p_norm - 1.0).abs() < DBL_EPSILON {
            // If p is already very close to a unit vector, avoid normalization
            p
        } else {
            p.normalize()
        };

        println!("Projection point before normalization: {:?}", p);
        println!("Projection point after normalization: {:?}", normalized_p);

        // For very small distances, we need to be careful with floating point precision
        // First check if the points are extremely close
        let direct_dist = x.0.sub(normalized_p).norm();
        if direct_dist < DBL_EPSILON {
            println!("Points are extremely close, returning minimal distance");
            return Angle(0.0);
        }

        let dist = x.distance(&Point(normalized_p));
        println!("Distance to projection: {}", dist.0);
        return dist;
    }

    // Returns a random orthonormal frame (three orthogonal unit-length vectors).
    fn random_frame() -> [Point; 3] {
        use rand::Rng;
        let mut rng = rand::thread_rng();

        // Generate a random point on the unit sphere.
        let dir = loop {
            let x = rng.gen_range(-1.0..1.0);
            let y = rng.gen_range(-1.0..1.0);
            let z = rng.gen_range(-1.0..1.0);
            let v = Vector::new(x, y, z);
            let norm2 = v.norm2();
            if norm2 > 0.0 && norm2 <= 1.0 {
                break Point(v.normalize());
            }
        };

        // Now choose a random orthogonal point.
        let ortho = dir.ortho();

        // Create a right-handed coordinate system.
        let third = Point(dir.0.cross(&ortho.0));

        let frame = [dir, ortho, third];

        println!("Random frame:");
        println!("  f[0] (p): {:?}, norm: {}", frame[0], frame[0].0.norm());
        println!("  f[1] (d1): {:?}, norm: {}", frame[1], frame[1].0.norm());
        println!("  f[2] (d2): {:?}, norm: {}", frame[2], frame[2].0.norm());
        println!("  Orthogonality check:");
        println!("    f[0]·f[1]: {}", frame[0].0.dot(&frame[1].0));
        println!("    f[0]·f[2]: {}", frame[0].0.dot(&frame[2].0));
        println!("    f[1]·f[2]: {}", frame[1].0.dot(&frame[2].0));

        frame
    }

    // Returns a random number in [0,1).
    fn random_float64() -> f64 {
        use rand::Rng;
        rand::thread_rng().gen()
    }

    // Returns true with 1/n probability.
    fn one_in(n: u32) -> bool {
        use rand::Rng;
        rand::thread_rng().gen_range(0..n) == 0
    }

    // Returns the maximum of two angles.
    fn max_angle(a: Angle, b: Angle) -> Angle {
        if a > b {
            a
        } else {
            b
        }
    }

    #[test]
    fn test_angle_contains_vertex() {
        let a = point_from_coords(1.0, 0.0, 0.0);
        let b = point_from_coords(0.0, 1.0, 0.0);
        let ref_b = b.reference_dir();

        // Degenerate angle ABA.
        assert!(
            !angle_contains_vertex(&a, &b, &a),
            "AngleContainsVertex({:?}, {:?}, {:?}) = true, want false",
            a,
            b,
            a
        );

        // An angle where A == referenceDir(B).
        assert!(
            angle_contains_vertex(&ref_b, &b, &a),
            "AngleContainsVertex({:?}, {:?}, {:?}) = false, want true",
            ref_b,
            b,
            a
        );

        // An angle where C == referenceDir(B).
        assert!(
            !angle_contains_vertex(&a, &b, &ref_b),
            "AngleContainsVertex({:?}, {:?}, {:?}) = true, want false",
            a,
            b,
            ref_b
        );

        // TODO: Add test for polygon tiling around a vertex once we have the RegularLoop implementation
    }

    #[test]
    fn test_edgeutil_intersection_error() {
        // We repeatedly construct two edges that cross near a random point "p", and
        // measure the distance from the actual intersection point "x" to the
        // exact intersection point and also to the edges.
        let distance_abs_error = Angle(3.0 * DBL_EPSILON);

        let mut max_point_dist = Angle(0.0);
        let mut max_edge_dist = Angle(0.0);

        // Reduce iterations for faster tests
        let iterations = if cfg!(debug_assertions) { 100 } else { 5000 };

        for _ in 0..iterations {
            // We construct two edges AB and CD that intersect near "p". The angle
            // between AB and CD (expressed as a slope) is chosen randomly between
            // 1e-15 and 1e15 such that its logarithm is uniformly distributed.
            // Similarly, two edge lengths approximately between 1e-15 and 1 are
            // chosen. The edge endpoints are chosen such that they are often very
            // close to the other edge (i.e., barely crossing). Taken together this
            // ensures that we test both long and very short edges that intersect at
            // both large and very small angles.
            //
            // Sometimes the edges we generate will not actually cross, in which case
            // we simply try again.
            let f = random_frame();
            let p = f[0];
            let d1 = f[1];
            let mut d2 = f[2];

            let slope = 1e-15 * (1e30_f64).powf(random_float64());
            println!("Slope used: {}", slope);
            d2 = Point(d1.0.add(d2.0.mul(slope)).normalize());

            // Find a pair of segments that cross.
            let (a, b, c, d) = loop {
                let ab_len = (1e-15_f64).powf(random_float64());
                let cd_len = (1e-15_f64).powf(random_float64());

                let mut a_fraction = (1e-5_f64).powf(random_float64());
                if one_in(2) {
                    a_fraction = 1.0 - a_fraction;
                }

                let mut c_fraction = (1e-5_f64).powf(random_float64());
                if one_in(2) {
                    c_fraction = 1.0 - c_fraction;
                }

                println!("ab_len: {}, cd_len: {}", ab_len, cd_len);
                println!("a_fraction: {}, c_fraction: {}", a_fraction, c_fraction);

                let a = Point(p.0.sub(d1.0.mul(a_fraction * ab_len)).normalize());
                let b = Point(p.0.add(d1.0.mul((1.0 - a_fraction) * ab_len)).normalize());
                let c = Point(p.0.sub(d2.0.mul(c_fraction * cd_len)).normalize());
                let d = Point(p.0.add(d2.0.mul((1.0 - c_fraction) * cd_len)).normalize());

                let mut crosser = EdgeCrosser::new(&a, &b);
                if crosser.crossing_sign(&c, &d) == Crossing::Cross {
                    println!("Found crossing segments:");
                    println!("  a: {:?}", a);
                    println!("  b: {:?}", b);
                    println!("  c: {:?}", c);
                    println!("  d: {:?}", d);
                    break (a, b, c, d);
                }
            };

            // Each constructed edge should be at most 1.5 * dblEpsilon away from the
            // original point P.
            let dist_ab = distance_from_segment(&p, &a, &b);
            println!("Distance calculation details:");
            println!("  p to a distance: {}", p.distance(&a).0);
            println!("  p to b distance: {}", p.distance(&b).0);
            println!("  |a-b|: {}", a.0.sub(b.0).norm());
            println!("  Computed dist_ab: {}", dist_ab.0);

            let want_dist = Angle(1.5 * DBL_EPSILON) + distance_abs_error;
            println!("  Want dist: {}", want_dist.0);

            assert!(
                dist_ab <= want_dist,
                "DistanceFromSegment({:?}, {:?}, {:?}) = {:.32}, want <= {:.32}",
                p,
                a,
                b,
                dist_ab.0,
                want_dist.0
            );

            let dist_cd = distance_from_segment(&p, &c, &d);
            assert!(
                dist_cd <= want_dist,
                "DistanceFromSegment({:?}, {:?}, {:?}) = {:.32}, want <= {:.32}",
                p,
                c,
                d,
                dist_cd.0,
                want_dist.0
            );

            // Verify that the expected intersection point is close to both edges and
            // also close to the original point P. (It might not be very close to P
            // if the angle between the edges is very small.)
            let expected = test_intersection_exact(&a, &b, &c, &d);

            let dist_expected_ab = distance_from_segment(&expected, &a, &b);
            let want_dist_expected = Angle(3.0 * DBL_EPSILON) + distance_abs_error;
            assert!(
                dist_expected_ab <= want_dist_expected,
                "DistanceFromSegment({:?}, {:?}, {:?}) = {:.32}, want <= {:.32}",
                expected,
                a,
                b,
                dist_expected_ab.0,
                want_dist_expected.0
            );

            let dist_expected_cd = distance_from_segment(&expected, &c, &d);
            assert!(
                dist_expected_cd <= want_dist_expected,
                "DistanceFromSegment({:?}, {:?}, {:?}) = {:.32}, want <= {:.32}",
                expected,
                c,
                d,
                dist_expected_cd.0,
                want_dist_expected.0
            );

            let dist_expected_p = expected.distance(&p);
            let want_dist_p = Angle(3.0 * DBL_EPSILON / slope) + Angle(INTERSECTION_ERROR);
            assert!(
                dist_expected_p <= want_dist_p,
                "{:?}.Distance({:?}) = {:.32}, want <= {:.32}",
                expected,
                p,
                dist_expected_p.0,
                want_dist_p.0
            );

            // Now we actually test the Intersection() method.
            let actual = intersection(&a, &b, &c, &d);
            let dist_ab = distance_from_segment(&actual, &a, &b);
            let dist_cd = distance_from_segment(&actual, &c, &d);
            let point_dist = expected.distance(&actual);

            let want_dist_intersection = Angle(INTERSECTION_ERROR) + distance_abs_error;
            assert!(
                dist_ab <= want_dist_intersection,
                "DistanceFromSegment({:?}, {:?}, {:?}) = {:.32}, want <= {:.32}",
                actual,
                a,
                b,
                dist_ab.0,
                want_dist_intersection.0
            );

            assert!(
                dist_cd <= want_dist_intersection,
                "DistanceFromSegment({:?}, {:?}, {:?}) = {:.32}, want <= {:.32}",
                actual,
                c,
                d,
                dist_cd.0,
                want_dist_intersection.0
            );

            let want_point_dist = Angle(INTERSECTION_ERROR);
            assert!(
                point_dist <= want_point_dist,
                "{:?}.Distance({:?}) = {:.32}, want <= {:.32}",
                expected,
                actual,
                point_dist.0,
                want_point_dist.0
            );

            max_edge_dist = max_angle(max_edge_dist, max_angle(dist_ab, dist_cd));
            max_point_dist = max_angle(max_point_dist, point_dist);
        }
    }
}
