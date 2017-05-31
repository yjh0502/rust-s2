/*
Copyright 2016 Google Inc. All rights reserved.
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

//! This file contains various predicates that are guaranteed to produce
//! correct, consistent results. They are also relatively efficient. This is
//! achieved by computing conservative error bounds and falling back to high
//! precision or even exact arithmetic when the result is uncertain. Such
//! predicates are useful in implementing robust algorithms.
//!
//! See also EdgeCrosser, which implements various exact
//! edge-crossing predicates more efficiently than can be done here.

use consts::*;
use s2::point::Point;

/// MAX_DETERMINANT_ERROR is the maximum error in computing (AxB).C where all vectors
/// are unit length. Using standard inequalities, it can be shown that
///
///  fl(AxB) = AxB + D where |D| <= (|AxB| + (2/sqrt(3))*|A|*|B|) * e
///
/// where "fl()" denotes a calculation done in floating-point arithmetic,
/// |x| denotes either absolute value or the L2-norm as appropriate, and
/// e is a reasonably small value near the noise level of floating point
/// number accuracy. Similarly,
///
///  fl(B.C) = B.C + d where |d| <= (|B.C| + 2*|B|*|C|) * e .
///
/// Applying these bounds to the unit-length vectors A,B,C and neglecting
/// relative error (which does not affect the sign of the result), we get
///
///  fl((AxB).C) = (AxB).C + d where |d| <= (3 + 2/sqrt(3)) * e
const MAX_DETERMINANT_ERROR: f64 = 1.8274 * DBL_EPSILON;

/// DET_ERROR_MULTIPLIER is the factor to scale the magnitudes by when checking
/// for the sign of set of points with certainty. Using a similar technique to
/// the one used for MAX_DETERMINANT_ERROR, the error is at most:
///
///   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
///
/// If the determinant magnitude is larger than this value then we know
/// its sign with certainty.
const DET_ERROR_MULTIPLIER: f64 = 3.2321 * DBL_EPSILON;

#[derive(PartialEq,Eq)]
pub enum Direction {
    Clockwise,
    Indeterminate,
    CounterClockwise,
}

/// sign returns true if the points A, B, C are strictly counterclockwise,
/// and returns false if the points are clockwise or collinear (i.e. if they are all
/// contained on some great circle).
///
/// Due to numerical errors, situations may arise that are mathematically
/// impossible, e.g. ABC may be considered strictly CCW while BCA is not.
/// However, the implementation guarantees the following:
///
/// If Sign(a,b,c), then !Sign(c,b,a) for all a,b,c.
pub fn sign(a: &Point, b: &Point, c: &Point) -> bool {
    // NOTE(dnadasi): In the C++ API the equivalent method here was known as "SimpleSign".

    // We compute the signed volume of the parallelepiped ABC. The usual
    // formula for this is (A ⨯ B) · C, but we compute it here using (C ⨯ A) · B
    // in order to ensure that ABC and CBA are not both CCW. This follows
    // from the following identities (which are true numerically, not just
    // mathematically):
    //
    //     (1) x ⨯ y == -(y ⨯ x)
    //     (2) -x · y == -(x · y)
    return c.0.cross(&a.0).dot(&b.0) > 0.;
}

/// robust_sign returns a Direction representing the ordering of the points.
/// CounterClockwise is returned if the points are in counter-clockwise order,
/// Clockwise for clockwise, and Indeterminate if any two points are the same (collinear),
/// or the sign could not completely be determined.
///
/// This function has additional logic to make sure that the above properties hold even
/// when the three points are coplanar, and to deal with the limitations of
/// floating-point arithmetic.
///
/// RobustSign satisfies the following conditions:
///
///  (1) RobustSign(a,b,c) == Indeterminate if and only if a == b, b == c, or c == a
///  (2) RobustSign(b,c,a) == RobustSign(a,b,c) for all a,b,c
///  (3) RobustSign(c,b,a) == -RobustSign(a,b,c) for all a,b,c
///
/// In other words:
///
///  (1) The result is Indeterminate if and only if two points are the same.
///  (2) Rotating the order of the arguments does not affect the result.
///  (3) Exchanging any two arguments inverts the result.
///
/// On the other hand, note that it is not true in general that
/// RobustSign(-a,b,c) == -RobustSign(a,b,c), or any similar identities
/// involving antipodal points.
pub fn robust_sign(a: &Point, b: &Point, c: &Point) -> Direction {
    let sign = triage_sign(a, b, c);
    if sign == Direction::Indeterminate {
        expensive_sign(a, b, c)
    } else {
        sign
    }
}


/// stable_sign reports the direction sign of the points in a numerically stable way.
/// Unlike triageSign, this method can usually compute the correct determinant sign
/// even when all three points are as collinear as possible. For example if three
/// points are spaced 1km apart along a random line on the Earth's surface using
/// the nearest representable points, there is only a 0.4% chance that this method
/// will not be able to find the determinant sign. The probability of failure
/// decreases as the points get closer together; if the collinear points are 1 meter
/// apart, the failure rate drops to 0.0004%.
///
/// This method could be extended to also handle nearly-antipodal points, but antipodal
/// points are rare in practice so it seems better to simply fall back to
/// exact arithmetic in that case.
pub fn stable_sign(a: &Point, b: &Point, c: &Point) -> Direction {
    let ab = &b.0 - &a.0;
    let ab2 = ab.norm2();
    let bc = &c.0 - &b.0;
    let bc2 = bc.norm2();
    let ca = &a.0 - &c.0;
    let ca2 = ca.norm2();

    // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been
    // cyclically permuted if necessary so that AB is the longest edge. (This
    // minimizes the magnitude of cross product.)  At the same time we also
    // compute the maximum error in the determinant.

    // The two shortest edges, pointing away from their common point.
    let (e1, e2, op) = if ab2 >= bc2 && ab2 >= ca2 {
        // AB is the longest edge.
        (ca, bc, &c.0)
    } else if bc2 >= ca2 {
        // BC is the longest edge.
        (ab, ca, &a.0)
    } else {
        // CA is the longest edge.
        (bc, ab, &b.0)
    };

    let det = -1. * e1.cross(&e2).dot(&op);
    let max_err = DET_ERROR_MULTIPLIER * (e1.norm2() * e2.norm2()).sqrt();

    // If the determinant isn't zero, within maxErr, we know definitively the point ordering.
    if det > max_err {
        Direction::CounterClockwise
    } else if det < -max_err {
        Direction::Clockwise
    } else {
        Direction::Indeterminate
    }
}

/// triage_sign returns the direction sign of the points. It returns Indeterminate if two
/// points are identical or the result is uncertain. Uncertain cases can be resolved, if
/// desired, by calling expensiveSign.
///
/// The purpose of this method is to allow additional cheap tests to be done without
/// calling expensiveSign.
pub fn triage_sign(a: &Point, b: &Point, c: &Point) -> Direction {
    let det = a.0.cross(&b.0).dot(&c.0);
    if det > MAX_DETERMINANT_ERROR {
        Direction::CounterClockwise
    } else if det < -MAX_DETERMINANT_ERROR {
        Direction::Clockwise
    } else {
        Direction::Indeterminate
    }
}

/// expensive_sign reports the direction sign of the points. It returns Indeterminate
/// if two of the input points are the same. It uses multiple-precision arithmetic
/// to ensure that its results are always self-consistent.
fn expensive_sign(a: &Point, b: &Point, c: &Point) -> Direction {
    // Return Indeterminate if and only if two points are the same.
    // This ensures RobustSign(a,b,c) == Indeterminate if and only if a == b, b == c, or c == a.
    // ie. Property 1 of RobustSign.
    if a == b || b == c || c == a {
        return Direction::Indeterminate;
    }

    // Next we try recomputing the determinant still using floating-point
    // arithmetic but in a more precise way. This is more expensive than the
    // simple calculation done by triageSign, but it is still *much* cheaper
    // than using arbitrary-precision arithmetic. This optimization is able to
    // compute the correct determinant sign in virtually all cases except when
    // the three points are truly collinear (e.g., three points on the equator).
    let det_sign = stable_sign(a, b, c);
    if det_sign != Direction::Indeterminate {
        det_sign
    } else {
        // Otherwise fall back to exact arithmetic and symbolic permutations.
        exact_sign(a, b, c, false)
    }
}

/// exact-sign reports the direction sign of the points using exact precision arithmetic.
fn exact_sign(_: &Point, _: &Point, _: &Point, _: bool) -> Direction {
    // In the C++ version, the final computation is performed using OpenSSL's
    // Bignum exact precision math library. The existence of an equivalent
    // library in Go is indeterminate. In C++, using the exact precision library
    // to solve this stage is ~300x slower than the above checks.
    // TODO(roberts): Select and incorporate an appropriate Go exact precision
    // floating point library for the remaining calculations.
    Direction::Indeterminate
}
