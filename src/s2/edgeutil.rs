/*
Copyright 2015 Google Inc. All rights reserved.

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

/*
package s2

import (
    "math"

    "github.com/golang/geo/r1"
    "github.com/golang/geo/r2"
    "github.com/golang/geo/r3"
    "github.com/golang/geo/s1"
)

const (
    // edgeClipErrorUVCoord is the maximum error in a u- or v-coordinate
    // compared to the exact result, assuming that the points A and B are in
    // the rectangle [-1,1]x[1,1] or slightly outside it (by 1e-10 or less).
    edgeClipErrorUVCoord = 2.25 * dblEpsilon

    // edgeClipErrorUVDist is the maximum distance from a clipped point to
    // the corresponding exact result. It is equal to the error in a single
    // coordinate because at most one coordinate is subject to error.
    edgeClipErrorUVDist = 2.25 * dblEpsilon

    // faceClipErrorRadians is the maximum angle between a returned vertex
    // and the nearest point on the exact edge AB. It is equal to the
    // maximum directional error in PointCross, plus the error when
    // projecting points onto a cube face.
    faceClipErrorRadians = 3 * dblEpsilon

    // faceClipErrorDist is the same angle expressed as a maximum distance
    // in (u,v)-space. In other words, a returned vertex is at most this far
    // from the exact edge AB projected into (u,v)-space.
    faceClipErrorUVDist = 9 * dblEpsilon

    // faceClipErrorUVCoord is the maximum angle between a returned vertex
    // and the nearest point on the exact edge AB expressed as the maximum error
    // in an individual u- or v-coordinate. In other words, for each
    // returned vertex there is a point on the exact edge AB whose u- and
    // v-coordinates differ from the vertex by at most this amount.
    faceClipErrorUVCoord = 9.0 * (1.0 / math.Sqrt2) * dblEpsilon

    // intersectsRectErrorUVDist is the maximum error when computing if a point
    // intersects with a given Rect. If some point of AB is inside the
    // rectangle by at least this distance, the result is guaranteed to be true;
    // if all points of AB are outside the rectangle by at least this distance,
    // the result is guaranteed to be false. This bound assumes that rect is
    // a subset of the rectangle [-1,1]x[-1,1] or extends slightly outside it
    // (e.g., by 1e-10 or less).
    intersectsRectErrorUVDist = 3 * math.Sqrt2 * dblEpsilon

    // intersectionError can be set somewhat arbitrarily, because the algorithm
    // uses more precision if necessary in order to achieve the specified error.
    // The only strict requirement is that intersectionError >= dblEpsilon
    // radians. However, using a larger error tolerance makes the algorithm more
    // efficient because it reduces the number of cases where exact arithmetic is
    // needed.
    intersectionError = s1.Angle(4 * dblEpsilon)

    // intersectionMergeRadius is used to ensure that intersection points that
    // are supposed to be coincident are merged back together into a single
    // vertex. This is required in order for various polygon operations (union,
    // intersection, etc) to work correctly. It is twice the intersection error
    // because two coincident intersection points might have errors in
    // opposite directions.
    intersectionMergeRadius = 2 * intersectionError
)
*/

use crate::consts::DBL_EPSILON;
use crate::s1::*;
use crate::s2::point::Point;
use crate::s2::predicates::sign;
use std::ops::{Mul, Sub};

// SimpleCrossing reports whether edge AB crosses CD at a point that is interior
// to both edges. Properties:
//
//  (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d)
//  (2) SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
pub fn simple_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> bool {
    // We compute the equivalent of Sign for triangles ACB, CBD, BDA,
    // and DAC. All of these triangles need to have the same orientation
    // (CW or CCW) for an intersection to exist.

    let ab = a.0.cross(&b.0);
    let acb = -(ab.dot(&c.0));
    let bda = ab.dot(&d.0);
    if acb * bda <= 0. {
        return false;
    }

    let cd = c.0.cross(&d.0);
    let cbd = -(cd.dot(&b.0));
    let dac = cd.dot(&a.0);

    (acb * cbd > 0.) && (acb * dac > 0.)
}

/// interpolate returns the point X along the line segment AB whose distance from A
/// is the given fraction "t" of the distance AB. Does NOT require that "t" be
/// between 0 and 1. Note that all distances are measured on the surface of
/// the sphere, so this is more complicated than just computing (1-t)*a + t*b
/// and normalizing the result.
pub fn interpolate(t: f64, a: &Point, b: &Point) -> Point {
    if t == 0f64 {
        return *a;
    }
    if t == 1f64 {
        return *b;
    }
    let ab = a.0.angle(&b.0).rad();
    interpolate_at_distance(&Angle::from(Rad(t * ab)), &a, &b)
}

/// interpolate_at_distance returns the point X along the line segment AB whose
/// distance from A is the angle ax.
pub fn interpolate_at_distance(ax: &Angle, a: &Point, b: &Point) -> Point {
    // aRad := ax.Radians()

    // Use cross to compute the tangent vector at A towards B. The
    // result is always perpendicular to A, even if A=B or A=-B, but it is not
    // necessarily unit length. (We effectively normalize it below.)
    let normal = a.cross(&b);
    let tangent = normal.0.cross(&a.0);

    // Now compute the appropriate linear combination of A and "tangent". With
    // infinite precision the result would always be unit length, but we
    // normalize it anyway to ensure that the error is within acceptable bounds.
    // (Otherwise errors can build up when the result of one interpolation is
    // fed into another interpolation.)
    let v = ((a.0 * ax.rad().cos()) + (tangent * (ax.rad().sin() / tangent.norm()))).normalize();
    Point(v)
}

// project returns the point along the edge AB that is closest to the point X.
// The fractional distance of this point along the edge AB can be obtained
// using distance_fraction.
//
// This requires that all points are unit length.
fn project(x: &Point, a: &Point, b: &Point) -> Point {
    let axb = &a.cross(b);
    // Find the closest point to X along the great circle through AB.
    let p = &x.sub(&axb.mul(x.0.dot(&axb.0) / axb.0.norm2()));

    // If this point is on the edge AB, then it's the closest point.
    if sign(axb, a, p) && sign(p, b, axb) {
        return p.normalize();
    }

    // Otherwise, the closest point is either A or B.
    if x.sub(a).0.norm2() <= x.sub(b).0.norm2() {
        return *a;
    }
    *b
}

// update_min_distance computes the distance from a point X to a line segment AB,
// and if either the distance was less than the given min_dist, or alwaysUpdate is
// true, the value and whether it was updated are returned.
fn update_min_distance(
    x: &Point,
    a: &Point,
    b: &Point,
    min_dist: ChordAngle,
    always_update: bool,
) -> (ChordAngle, bool) {
    let (d, ok) = interior_dist(x, a, b, min_dist, always_update);
    if ok {
        // Minimum distance is attained along the edge interior.
        return (d, true);
    }
    // Otherwise the minimum distance is to one of the endpoints.
    let (xa2, xb2) = (x.sub(a).0.norm2(), x.sub(b).0.norm2());
    let dist = ChordAngle(xa2.min(xb2));
    if !always_update && dist >= min_dist {
        return (min_dist, false);
    }
    (dist, true)
}

// update_max_distance checks if the distance from X to the edge AB is greater
// than maxDist, and if so, returns the updated value and true.
// Otherwise it returns false. The case A == B is handled correctly.
fn update_max_distance(
    x: &Point,
    a: &Point,
    b: &Point,
    max_dist: ChordAngle,
) -> (ChordAngle, bool) {
    let mut dist = Point::chordangle(x, a).max(Point::chordangle(x, b));

    if dist > chordangle::RIGHT {
        let (dist2, _) = update_min_distance(&Point(x.0.mul(-1.)), a, b, dist, true);
        dist = chordangle::STRAIGHT - dist2;
    }
    if max_dist < dist {
        return (dist, true);
    }
    (max_dist, false)
}

// distance_from_segment returns the distance of point X from line segment AB.
// The points are expected to be normalized. The result is very accurate for small
// distances but may have some numerical error if the distance is large
// (approximately pi/2 or greater). The case A == B is handled correctly.
pub fn distance_from_segment(x: &Point, a: &Point, b: &Point) -> Angle {
    let (min_dist, _) = update_min_distance(x, a, b, ChordAngle::default(), true);
    Angle::from(min_dist)
}

// interior_dist returns the shortest distance from point x to edge ab, assuming
// that the closest point to X is interior to AB. If the closest point is not
// interior to AB, interiorDist returns (minDist, false). If alwaysUpdate is set to
// false, the distance is only updated when the value exceeds certain the given minDist.
fn interior_dist(
    x: &Point,
    a: &Point,
    b: &Point,
    min_dist: ChordAngle,
    always_update: bool,
) -> (ChordAngle, bool) {
    // Chord distance of x to both end points a and b.
    let (xa2, xb2) = ((x.sub(a)).0.norm2(), x.sub(b).0.norm2());
    // The closest point on AB could either be one of the two vertices (the
    // vertex case) or in the interior (the interior case). Let C = A x B.
    // If X is in the spherical wedge extending from A to B around the axis
    // through C, then we are in the interior case. Otherwise we are in the
    // vertex case.
    //
    // Check whether we might be in the interior case. For this to be true, XAB
    // and XBA must both be acute angles. Checking this condition exactly is
    // expensive, so instead we consider the planar triangle ABX (which passes
    // through the sphere's interior). The planar angles XAB and XBA are always
    // less than the corresponding spherical angles, so if we are in the
    // interior case then both of these angles must be acute.
    //
    // We check this by computing the squared edge lengths of the planar
    // triangle ABX, and testing whether angles XAB and XBA are both acute using
    // the law of cosines:
    //
    //            | XA^2 - XB^2 | < AB^2      (*)
    //
    // This test must be done conservatively (taking numerical errors into
    // account) since otherwise we might miss a situation where the true minimum
    // distance is achieved by a point on the edge interior.
    //
    // There are two sources of error in the expression above (*).  The first is
    // that points are not normalized exactly; they are only guaranteed to be
    // within 2 * DBL_EPSILON of unit length.  Under the assumption that the two
    // sides of (*) are nearly equal, the total error due to normalization errors
    // can be shown to be at most
    //
    //        2 * DBL_EPSILON * (XA^2 + XB^2 + AB^2) + 8 * DBL_EPSILON ^ 2 .
    //
    // The other source of error is rounding of results in the calculation of (*).
    // Each of XA^2, XB^2, AB^2 has a maximum relative error of 2.5 * DBL_EPSILON,
    // plus an additional relative error of 0.5 * DBL_EPSILON in the final
    // subtraction which we further bound as 0.25 * DBL_EPSILON * (XA^2 + XB^2 +
    // AB^2) for convenience.  This yields a final error bound of
    //
    //        4.75 * DBL_EPSILON * (XA^2 + XB^2 + AB^2) + 8 * DBL_EPSILON ^ 2 .
    let ab2 = a.sub(b).0.norm2();
    let max_error = 4.75 * DBL_EPSILON * (xa2 + xb2 + ab2) + 8. * DBL_EPSILON * DBL_EPSILON;
    if (xa2 - xb2).abs() >= ab2 + max_error {
        return (min_dist, false);
    }

    // The minimum distance might be to a point on the edge interior. Let R
    // be closest point to X that lies on the great circle through AB. Rather
    // than computing the geodesic distance along the surface of the sphere,
    // instead we compute the "chord length" through the sphere's interior.
    //
    // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q
    // is the point X projected onto the plane through the great circle AB.
    // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
    // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it
    // is faster and the corresponding distance on the Earth's surface is
    // accurate to within 1% for distances up to about 1800km.
    let c = &a.cross(b).0;
    let c2 = c.norm2();
    let x_dot_c = x.0.dot(c);
    let x_dot_c2 = x_dot_c * x_dot_c;
    if !always_update && x_dot_c2 > c2 * min_dist.0 {
        // The closest point on the great circle AB is too far away.  We need to
        // test this using ">" rather than ">=" because the actual minimum bound
        // on the distance is (xDotC2 / c2), which can be rounded differently
        // than the (more efficient) multiplicative test above.
        return (min_dist, false);
    }

    // Otherwise we do the exact, more expensive test for the interior case.
    // This test is very likely to succeed because of the conservative planar
    // test we did initially.
    //
    // TODO(roberts): Ensure that the errors in test are accurately reflected in the
    // minUpdateInteriorDistanceMaxError.
    let cx = &c.cross(&x.0);
    if a.sub(x).0.dot(cx) >= 0. || b.sub(x).0.dot(cx) <= 0. {
        return (min_dist, false);
    }

    // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
    // This calculation has good accuracy for all chord lengths since it
    // is based on both the dot product and cross product (rather than
    // deriving one from the other). However, note that the chord length
    // representation itself loses accuracy as the angle approaches π.
    let qr = 1. - (cx.norm2() / c2).sqrt();
    let dist = ChordAngle((x_dot_c2 / c2) + (qr * qr));

    if !always_update && dist >= min_dist {
        return (min_dist, false);
    }
    return (dist, true);
}

#[cfg(test)]
#[allow(non_upper_case_globals)]
mod tests {
    use super::*;
    use std::f64::consts::{PI, SQRT_2};

    use crate::r3::vector::Vector;

    fn f64_near(x: f64, y: f64, epsilon: f64) -> bool {
        (x - y).abs() <= epsilon
    }

    #[test]
    fn test_edge_distances_check_distance() {
        let tests = [
            (
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                0.,
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                0.,
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 1.,
                    y: 3.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                0.,
                Vector {
                    x: 1.,
                    y: 3.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 0.,
                    y: 0.,
                    z: 1.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                PI / 2.,
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 0.,
                    y: 0.,
                    z: -1.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                PI / 2.,
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: -1.,
                    y: -1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                0.75 * PI,
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                PI / 4.,
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 0.,
                    y: -1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                PI / 2.,
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 0.,
                    y: -1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: -1.,
                    y: 1.,
                    z: 0.,
                },
                PI / 2.,
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: -1.,
                    y: -1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: -1.,
                    y: 1.,
                    z: 0.,
                },
                PI / 2.,
                Vector {
                    x: -1.,
                    y: 1.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 1.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                f64::sqrt(1.0 / 3.0).asin(),
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 1.,
                    y: 1.,
                    z: -1.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                f64::sqrt(1.0 / 3.0).asin(),
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: -1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                0.75 * PI,
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: 0.,
                    y: 0.,
                    z: -1.,
                },
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                PI / 2.,
                Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
            ),
            (
                Vector {
                    x: -1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                PI,
                Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
            ),
        ];

        for &(x, a, b, dist_rad, want) in &tests {
            let x = Point(x.normalize());
            let a = Point(a.normalize());
            let b = Point(b.normalize());
            let want_p = Point(want.normalize());

            let d = distance_from_segment(&x, &a, &b).rad();
            assert!(f64_near(d, dist_rad, 1e-15));

            let closest = project(&x, &a, &b);
            assert!(closest.approx_eq(&want_p));

            let (_min_distance, ok) = update_min_distance(&x, &a, &b, ChordAngle(0.), false);
            assert!(!ok);

            let (min_distance, ok) = update_min_distance(&x, &a, &b, ChordAngle::inf(), false);
            assert!(ok);

            assert!(f64_near(dist_rad, Angle::from(min_distance).rad(), 1e-15));
        }
    }

    #[test]
    fn test_edge_distances_update_min_interior_distance_lower_bound_optimization_is_conservative() {
        // Verifies that alwaysUpdateMinInteriorDistance computes the lower bound
        // on the true distance conservatively.  (This test used to fail.)
        let x = Point::from_coords(
            -0.017952729194524016,
            -0.30232422079175203,
            0.95303607751077712,
        );
        let a = Point::from_coords(
            -0.017894725505830295,
            -0.30229974986194175,
            0.95304493075220664,
        );
        let b = Point::from_coords(
            -0.017986591360900289,
            -0.30233851195954353,
            0.95303090543659963,
        );

        let (min_distance, ok) = update_min_distance(&x, &a, &b, ChordAngle::inf(), false);
        assert!(ok);
        let min_distance = min_distance.successor();
        let (_, ok) = update_min_distance(&x, &a, &b, min_distance, false);
        assert!(ok);
    }

    #[test]
    fn test_edge_distances_update_min_interior_distance_rejection_test_is_conservative() {
        // This test checks several representative cases where previously
        // UpdateMinInteriorDistance was failing to update the distance because a
        // rejection test was not being done conservatively.
        //
        // Note that all of the edges AB in this test are nearly antipodal.

        let min_dist = ChordAngle::from_squared_length(6.3897233584120815e-26);

        let tests = [
            (
                Point(Vector {
                    x: 1.,
                    y: -4.6547732744037044e-11,
                    z: -5.6374428459823598e-89,
                }),
                Point(Vector {
                    x: 1.,
                    y: -8.9031850507928352e-11,
                    z: 0.,
                }),
                Point(Vector {
                    x: -0.99999999999996347,
                    y: 2.7030110029169596e-07,
                    z: 1.555092348806121e-99,
                }),
                min_dist,
                false,
            ),
            (
                Point(Vector {
                    x: 1.,
                    y: -4.7617930898495072e-13,
                    z: 0.,
                }),
                Point(Vector {
                    x: -1.,
                    y: -1.6065916409055676e-10,
                    z: 0.,
                }),
                Point(Vector {
                    x: 1.,
                    y: 0.,
                    z: 9.9964883247706732e-35,
                }),
                min_dist,
                false,
            ),
            (
                Point(Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -8.4965026896454536e-11,
                    z: 0.,
                }),
                Point(Vector {
                    x: -0.99999999999966138,
                    y: 8.2297529603339328e-07,
                    z: 9.6070344113320997e-21,
                }),
                min_dist,
                false,
            ),
        ];
        for (x, a, b, min_dist, _want) in &tests {
            let (_, ok) = update_min_distance(x, a, b, *min_dist, false);
            assert!(ok);
        }
    }

    #[test]
    fn test_edge_distances_check_max_distance() {
        struct Test {
            x: Vector,
            a: Vector,
            b: Vector,
            dist_rad: f64,
        };
        let tests: [Test; 12] = [
            Test {
                x: Vector {
                    x: 1.,
                    y: 0.,
                    z: 1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                dist_rad: PI / 2.,
            },
            Test {
                x: Vector {
                    x: 1.,
                    y: 0.,
                    z: -1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                dist_rad: PI / 2.,
            },
            Test {
                x: Vector {
                    x: 0.,
                    y: 1.,
                    z: 1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                dist_rad: PI / 2.,
            },
            Test {
                x: Vector {
                    x: 0.,
                    y: 1.,
                    z: -1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                dist_rad: PI / 2.,
            },
            Test {
                x: Vector {
                    x: 1.,
                    y: 1.,
                    z: 1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                dist_rad: f64::sqrt(2. / 3.).asin(),
            },
            Test {
                x: Vector {
                    x: 1.,
                    y: 1.,
                    z: -1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                dist_rad: f64::sqrt(2. / 3.).asin(),
            },
            Test {
                x: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                a: Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                b: Vector {
                    x: 1.,
                    y: -1.,
                    z: 0.,
                },
                dist_rad: PI / 4.,
            },
            Test {
                x: Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                },
                a: Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                b: Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                },
                dist_rad: PI / 4.,
            },
            Test {
                x: Vector {
                    x: 0.,
                    y: 0.,
                    z: 1.,
                },
                a: Vector {
                    x: 0.,
                    y: 1.,
                    z: 1.,
                },
                b: Vector {
                    x: 0.,
                    y: -1.,
                    z: 1.,
                },
                dist_rad: PI / 4.,
            },
            Test {
                x: Vector {
                    x: 0.,
                    y: 0.,
                    z: 1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 1.,
                    y: 0.,
                    z: -1.,
                },
                dist_rad: 3. * PI / 4.,
            },
            Test {
                x: Vector {
                    x: 0.,
                    y: 0.,
                    z: 1.,
                },
                a: Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                },
                b: Vector {
                    x: 1.,
                    y: 1.,
                    z: -SQRT_2,
                },
                dist_rad: 3. * PI / 4.,
            },
            Test {
                x: Vector {
                    x: 0.,
                    y: 0.,
                    z: 1.,
                },
                a: Vector {
                    x: 0.,
                    y: 0.,
                    z: -1.,
                },
                b: Vector {
                    x: 0.,
                    y: 0.,
                    z: -1.,
                },
                dist_rad: PI,
            },
        ];

        for test in &tests {
            let x = &Point(test.x.normalize());
            let a = &Point(test.a.normalize());
            let b = &Point(test.b.normalize());

            let max_distance = chordangle::STRAIGHT;
            let (_, ok) = update_max_distance(x, a, b, max_distance);
            assert!(!ok);

            let max_distance = chordangle::NEGATIVE;
            let (max_distance, ok) = update_max_distance(x, a, b, max_distance);
            assert!(ok);

            assert!(f64_near(
                test.dist_rad,
                Angle::from(max_distance).rad(),
                1e-15
            ));
        }
    }
}
