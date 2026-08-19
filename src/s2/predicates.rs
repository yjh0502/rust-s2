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
//!
//! ## Porting notes (functions and tests, source and verification status)
//!
//! This file was originally ported from Go's `s2/predicates.go`
//! (golang/geo). The notes below track, for each piece touched or added
//! while fixing https://github.com/yjh0502/rust-s2/issues/132, which
//! upstream(s) it was checked against and whether that check found the
//! two upstreams in agreement or not.
//!
//! - `exact_sign` / `symbolically_perturbed_sign`: ported from Go's
//!   `exactSign` / `symbolicallyPerturbedSign`, then independently checked
//!   line-by-line against C++'s `ExactSign` / `SymbolicallyPerturbedSign`
//!   (`s2predicates.cc`, using `ExactFloat`/Bignum where Go uses
//!   `big.Float` and we use `bigdecimal::BigDecimal`). **Verified
//!   identical** in both: same branch order, same cascade, no
//!   discrepancy.
//! - `stable_sign`: pre-existing (Go-derived) function; while adding the
//!   test below, checked against both C++'s `StableSign` and Go's
//!   `stableSign`. **Discrepancy found**: both upstreams guard against the
//!   max-error computation underflowing to zero (`minNoUnderflowError` /
//!   `kMinNoUnderflowError`) before trusting the result; this Rust port
//!   was missing that guard entirely. Fixed to match both (they agree
//!   with each other, so there was no C++-vs-Go conflict to resolve).
//! - `test_sign`, `test_robust_sign`, `test_robust_sign_equalities`:
//!   ported from Go's `TestPredicatesSign` / `TestPredicatesRobustSign` /
//!   `TestPredicatesRobustSignEqualities`, which existed in this file as
//!   dead, commented-out Go source before this fix. Not cross-checked
//!   against C++'s `s2predicates_test.cc` test bodies (which structure
//!   the same cases differently, as `Sign.CollinearPoints`), since the
//!   Go-derived tests already cover the same points and directions.
//! - `test_symbolic_perturbation_code_coverage`: ported from C++'s
//!   `Sign.SymbolicPerturbationCodeCoverage` (`s2predicates_test.cc`).
//!   No Go equivalent exists. All 13 cases pass, which is independent
//!   confirmation that `symbolically_perturbed_sign` above is correct.
//! - `test_stable_sign_underflow`: ported from C++'s
//!   `Sign.StableSignUnderflow` (`s2predicates_test.cc`). No Go
//!   equivalent exists. This is the regression test for the `stable_sign`
//!   discrepancy noted above.
//! - Not ported: C++'s `SignTest::StressTest`, a randomized great-circle
//!   fuzz test with substantial supporting fixture code; and Go's/C++'s
//!   `CompareDistances`, `CompareDistance`, `SignDotProd`,
//!   `CircleEdgeIntersectionOrdering`, which back closest/furthest-edge
//!   query functionality not yet ported to this crate (see
//!   https://github.com/yjh0502/rust-s2/issues/133).

use crate::consts::*;
use crate::r3::precisevector::PreciseVector;
use crate::s2::point::Point;
use bigdecimal::BigDecimal;
use bigdecimal::num_bigint::Sign as BigSign;

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

#[derive(PartialEq, Eq, Clone, Copy, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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
    c.0.cross(&a.0).dot(&b.0) > 0f64
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
    let ab = b.0 - a.0;
    let ab2 = ab.norm2();
    let bc = c.0 - b.0;
    let bc2 = bc.norm2();
    let ca = a.0 - c.0;
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

    let det = -e1.cross(&e2).dot(op);
    let max_err = DET_ERROR_MULTIPLIER * (e1.norm2() * e2.norm2()).sqrt();

    // Errors smaller than this value may not be accurate due to underflow.
    let min_no_underflow_error = DET_ERROR_MULTIPLIER * f64::MIN_POSITIVE.sqrt();
    if max_err < min_no_underflow_error {
        return Direction::Indeterminate;
    }

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
        exact_sign(a, b, c, true)
    }
}

/// direction_sign converts the sign of a high-precision value into a Direction,
/// following the convention that CounterClockwise corresponds to a positive
/// sign, Clockwise to a negative sign, and Indeterminate to zero.
fn direction_sign(v: &BigDecimal) -> Direction {
    match v.sign() {
        BigSign::Plus => Direction::CounterClockwise,
        BigSign::Minus => Direction::Clockwise,
        BigSign::NoSign => Direction::Indeterminate,
    }
}

/// flip reverses a Direction, leaving Indeterminate unchanged.
fn flip(d: Direction) -> Direction {
    match d {
        Direction::Clockwise => Direction::CounterClockwise,
        Direction::CounterClockwise => Direction::Clockwise,
        Direction::Indeterminate => Direction::Indeterminate,
    }
}

/// exact_sign reports the direction sign of the points using exact precision
/// arithmetic, falling back to symbolic perturbation if the three points are
/// (exactly) collinear and `perturb` is true.
fn exact_sign(a: &Point, b: &Point, c: &Point, perturb: bool) -> Direction {
    // Sort the three points in lexicographic order, keeping track of the
    // sign of the permutation. (Each exchange inverts the sign of the
    // determinant.)
    let mut perm_sign = Direction::CounterClockwise;
    let (mut pa, mut pb, mut pc) = (a, b, c);
    if pa.0 > pb.0 {
        std::mem::swap(&mut pa, &mut pb);
        perm_sign = flip(perm_sign);
    }
    if pb.0 > pc.0 {
        std::mem::swap(&mut pb, &mut pc);
        perm_sign = flip(perm_sign);
    }
    if pa.0 > pb.0 {
        std::mem::swap(&mut pa, &mut pb);
        perm_sign = flip(perm_sign);
    }

    // Construct multiple-precision versions of the sorted points and
    // compute their exact 3x3 determinant.
    let xa = PreciseVector::from(pa.0);
    let xb = PreciseVector::from(pb.0);
    let xc = PreciseVector::from(pc.0);
    let xb_cross_xc = xb.cross(xc.clone());
    let det = xa.dot(xb_cross_xc.clone());

    // Unlike C++'s Bignum or Go's big.Float, BigDecimal performs +, -, and *
    // exactly (no rounding), so the sign of `det` is the true sign of the
    // determinant and no precision analysis is needed here.

    // If the exact determinant is non-zero, we're done.
    let mut det_sign = direction_sign(&det);
    if det_sign == Direction::Indeterminate && perturb {
        // The points are exactly collinear. For example, this happens when
        // three points are spaced at exactly equal intervals along a
        // longitude line. Resolve the sign via symbolic perturbation. This
        // case does not happen in practice because IEEE 754 rounding errors
        // prevent three exactly collinear points from happening, other than
        // in tests that explicitly construct such situations.
        det_sign = symbolically_perturbed_sign(&xa, &xb, &xc, &xb_cross_xc);
    }

    // permSign is always CounterClockwise or Clockwise (never Indeterminate),
    // so the combined sign is Indeterminate exactly when det_sign is.
    if det_sign == Direction::Indeterminate {
        Direction::Indeterminate
    } else if perm_sign == det_sign {
        Direction::CounterClockwise
    } else {
        Direction::Clockwise
    }
}

/// symbolically_perturbed_sign reports the sign of the determinant of three
/// points A, B, C under a model where every possible Point is slightly
/// perturbed by a unique infinitesimal amount such that no three perturbed
/// points are collinear and no four are coplanar. The perturbations are so
/// small that they do not change the sign of any determinant that was
/// non-zero before the perturbation, and therefore can be safely ignored
/// unless the exact determinant of the three points is zero.
///
/// This returns CounterClockwise or Clockwise according to the sign of the
/// determinant after the symbolic perturbations are taken into account.
///
/// Since the symbolic perturbation of a given point is fixed (i.e. it does
/// not depend on the other two arguments), the results of this function are
/// always self-consistent: it will never return a result that corresponds to
/// an impossible configuration of non-degenerate points.
///
/// Requires that the exact 3x3 determinant of a, b, c is zero, and that the
/// points are distinct with a < b < c in lexicographic order.
///
/// Reference: "Simulation of Simplicity" (Edelsbrunner and Muecke, ACM
/// Transactions on Graphics, 1990).
fn symbolically_perturbed_sign(
    a: &PreciseVector,
    b: &PreciseVector,
    c: &PreciseVector,
    b_cross_c: &PreciseVector,
) -> Direction {
    // This is a direct translation of the fifteen-case expansion of the
    // fully-perturbed 3x3 determinant of A, B, C in powers of the
    // perturbation parameter; see the Go/C++ implementations for the
    // derivation.
    let mut d = direction_sign(b_cross_c.z()); // da.Z
    if d != Direction::Indeterminate {
        return d;
    }
    d = direction_sign(b_cross_c.y()); // da.Y
    if d != Direction::Indeterminate {
        return d;
    }
    d = direction_sign(b_cross_c.x()); // da.X
    if d != Direction::Indeterminate {
        return d;
    }

    d = direction_sign(&(c.x() * a.y() - c.y() * a.x())); // db.Z
    if d != Direction::Indeterminate {
        return d;
    }
    d = direction_sign(c.x()); // db.Z * da.Y
    if d != Direction::Indeterminate {
        return d;
    }
    d = flip(direction_sign(c.y())); // db.Z * da.X
    if d != Direction::Indeterminate {
        return d;
    }

    d = direction_sign(&(c.z() * a.x() - c.x() * a.z())); // db.Y
    if d != Direction::Indeterminate {
        return d;
    }
    d = direction_sign(c.z()); // db.Y * da.X
    if d != Direction::Indeterminate {
        return d;
    }

    d = direction_sign(&(a.x() * b.y() - a.y() * b.x())); // dc.Z
    if d != Direction::Indeterminate {
        return d;
    }
    d = flip(direction_sign(b.x())); // dc.Z * da.Y
    if d != Direction::Indeterminate {
        return d;
    }
    d = direction_sign(b.y()); // dc.Z * da.X
    if d != Direction::Indeterminate {
        return d;
    }
    d = direction_sign(a.x()); // dc.Z * db.Y
    if d != Direction::Indeterminate {
        return d;
    }
    Direction::CounterClockwise // dc.Z * db.Y * da.X
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r3::vector::Vector;

    fn pt(x: f64, y: f64, z: f64) -> Point {
        Point(Vector { x, y, z })
    }

    struct SignTest {
        p1: Point,
        p2: Point,
        p3: Point,
        want: bool,
    }

    #[test]
    fn test_sign() {
        let tests = vec![
            SignTest {
                p1: pt(1., 0., 0.),
                p2: pt(0., 1., 0.),
                p3: pt(0., 0., 1.),
                want: true,
            },
            SignTest {
                p1: pt(0., 1., 0.),
                p2: pt(0., 0., 1.),
                p3: pt(1., 0., 0.),
                want: true,
            },
            SignTest {
                p1: pt(0., 0., 1.),
                p2: pt(1., 0., 0.),
                p3: pt(0., 1., 0.),
                want: true,
            },
            SignTest {
                p1: pt(1., 1., 0.),
                p2: pt(0., 1., 1.),
                p3: pt(1., 0., 1.),
                want: true,
            },
            SignTest {
                p1: pt(-3., -1., 4.),
                p2: pt(2., -1., -3.),
                p3: pt(1., -2., 0.),
                want: true,
            },
            // All degenerate cases of sign(). Let M_1, M_2, ... be the sequence of
            // submatrices whose determinant sign is tested by that function. Then the
            // i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
            //
            //    det(M) = 0
            //    det(M_j) = 0 for j < i
            //    det(M_i) != 0
            //    A < B < C in lexicographic order.
            // det(M_1) = b0*c1 - b1*c0
            SignTest {
                p1: pt(-3., -1., 0.),
                p2: pt(-2., 1., 0.),
                p3: pt(1., -2., 0.),
                want: false,
            },
            // det(M_2) = b2*c0 - b0*c2
            SignTest {
                p1: pt(-6., 3., 3.),
                p2: pt(-4., 2., -1.),
                p3: pt(-2., 1., 4.),
                want: false,
            },
            // det(M_3) = b1*c2 - b2*c1
            SignTest {
                p1: pt(0., -1., -1.),
                p2: pt(0., 1., -2.),
                p3: pt(0., 2., 1.),
                want: false,
            },
            // From this point onward, B or C must be zero, or B is proportional to C.
            // det(M_4) = c0*a1 - c1*a0
            SignTest {
                p1: pt(-1., 2., 7.),
                p2: pt(2., 1., -4.),
                p3: pt(4., 2., -8.),
                want: false,
            },
            // det(M_5) = c0
            SignTest {
                p1: pt(-4., -2., 7.),
                p2: pt(2., 1., -4.),
                p3: pt(4., 2., -8.),
                want: false,
            },
            // det(M_6) = -c1
            SignTest {
                p1: pt(0., -5., 7.),
                p2: pt(0., -4., 8.),
                p3: pt(0., -2., 4.),
                want: false,
            },
            // det(M_7) = c2*a0 - c0*a2
            SignTest {
                p1: pt(-5., -2., 7.),
                p2: pt(0., 0., -2.),
                p3: pt(0., 0., -1.),
                want: false,
            },
            // det(M_8) = c2
            SignTest {
                p1: pt(0., -2., 7.),
                p2: pt(0., 0., 1.),
                p3: pt(0., 0., 2.),
                want: false,
            },
        ];

        for t in &tests {
            let got = sign(&t.p1, &t.p2, &t.p3);
            assert_eq!(
                got, t.want,
                "sign({:?}, {:?}, {:?}) = {}, want {}",
                t.p1, t.p2, t.p3, got, t.want
            );
            if t.want {
                // For these cases we can test the reversibility condition.
                let reversed = sign(&t.p3, &t.p2, &t.p1);
                assert_ne!(
                    reversed, t.want,
                    "sign({:?}, {:?}, {:?}) = {}, want {}",
                    t.p3, t.p2, t.p1, reversed, !t.want
                );
            }
        }
    }

    // The following points happen to be *exactly collinear* along a line that
    // is approximately tangent to the surface of the unit sphere. In fact, C
    // is the exact midpoint of the line segment AB. All of these points are
    // close enough to unit length to satisfy Vector::is_unit().
    fn po_a() -> Point {
        pt(
            0.72571927877036835,
            0.46058825605889098,
            0.51106749730504852,
        )
    }
    fn po_b() -> Point {
        pt(0.7257192746638208, 0.46058826573818168, 0.51106749441312738)
    }
    fn po_c() -> Point {
        pt(
            0.72571927671709457,
            0.46058826089853633,
            0.51106749585908795,
        )
    }

    // The points "x1" and "x2" are exactly proportional, i.e. they both lie
    // on a common line through the origin. Both points are considered to be
    // normalized, and in fact they both satisfy (x == x.normalize()).
    // Therefore the triangle (x1, x2, -x1) consists of three distinct points
    // that all lie on a common line through the origin.
    fn x1() -> Point {
        pt(0.99999999999999989, 1.4901161193847655e-08, 0.)
    }
    fn x2() -> Point {
        pt(1., 1.4901161193847656e-08, 0.)
    }

    // Here are two more points that are distinct, exactly proportional, and
    // that satisfy (x == x.normalize()).
    fn x3() -> Point {
        Point(
            Vector {
                x: 1.,
                y: 1.,
                z: 1.,
            }
            .normalize(),
        )
    }
    fn x4() -> Point {
        Point(x3().0 * 0.99999999999999989)
    }

    // The following three points demonstrate that normalize() is not
    // idempotent, i.e. y0.normalize() != y0.normalize().normalize(). Both
    // points are exactly proportional.
    fn y0() -> Point {
        pt(1., 1., 0.)
    }
    fn y1() -> Point {
        Point(y0().0.normalize())
    }
    fn y2() -> Point {
        Point(y1().0.normalize())
    }

    #[test]
    fn test_robust_sign_equalities() {
        let tests = vec![
            (Point(po_c().0 - po_a().0), Point(po_b().0 - po_c().0), true),
            (x1(), Point(x1().0.normalize()), true),
            (x2(), Point(x2().0.normalize()), true),
            (x3(), Point(x3().0.normalize()), true),
            (x4(), Point(x4().0.normalize()), true),
            (x3(), x4(), false),
            (y1(), y2(), false),
            (y2(), Point(y2().0.normalize()), true),
        ];

        for (p1, p2, want) in tests {
            assert_eq!(
                p1.0 == p2.0,
                want,
                "testing equality for robust_sign: {:?} == {:?}",
                p1,
                p2
            );
        }
    }

    #[test]
    fn test_robust_sign() {
        let x = pt(1., 0., 0.);
        let y = pt(0., 1., 0.);
        let z = pt(0., 0., 1.);

        struct RobustSignTest {
            p1: Point,
            p2: Point,
            p3: Point,
            want: Direction,
        }

        let tests = vec![
            // Simple collinear points test cases.
            // a == b != c
            RobustSignTest {
                p1: x,
                p2: x,
                p3: z,
                want: Direction::Indeterminate,
            },
            // a != b == c
            RobustSignTest {
                p1: x,
                p2: y,
                p3: y,
                want: Direction::Indeterminate,
            },
            // c == a != b
            RobustSignTest {
                p1: z,
                p2: x,
                p3: z,
                want: Direction::Indeterminate,
            },
            // CCW
            RobustSignTest {
                p1: x,
                p2: y,
                p3: z,
                want: Direction::CounterClockwise,
            },
            // CW
            RobustSignTest {
                p1: z,
                p2: y,
                p3: x,
                want: Direction::Clockwise,
            },
            // Edge cases: The following points happen to be *exactly collinear*
            // along a line that is approximately tangent to the surface of the
            // unit sphere. In fact, C is the exact midpoint of the line segment
            // AB. All of these points are close enough to unit length to satisfy
            // Vector::is_unit(). This used to only resolve to Indeterminate
            // before exact_sign() was implemented.
            RobustSignTest {
                p1: po_a(),
                p2: po_b(),
                p3: po_c(),
                want: Direction::Clockwise,
            },
            // The points "x1" and "x2" are exactly proportional, i.e. they both
            // lie on a common line through the origin. Both points are considered
            // to be normalized, and in fact they both satisfy (x ==
            // x.normalize()). Therefore the triangle (x1, x2, -x1) consists of
            // three distinct points that all lie on a common line through the
            // origin.
            RobustSignTest {
                p1: x1(),
                p2: x2(),
                p3: Point(x1().0 * -1.0),
                want: Direction::CounterClockwise,
            },
            // Here are two more points that are distinct, exactly proportional,
            // and that satisfy (x == x.normalize()).
            RobustSignTest {
                p1: x3(),
                p2: x4(),
                p3: Point(x3().0 * -1.0),
                want: Direction::Clockwise,
            },
            // The following points demonstrate that normalize() is not
            // idempotent, i.e. y0.normalize() != y0.normalize().normalize().
            // Both points satisfy is_unit(), though, and the two points are
            // exactly proportional.
            RobustSignTest {
                p1: y1(),
                p2: y2(),
                p3: Point(y1().0 * -1.0),
                want: Direction::CounterClockwise,
            },
        ];

        for t in &tests {
            let result = robust_sign(&t.p1, &t.p2, &t.p3);
            assert_eq!(
                result, t.want,
                "robust_sign({:?}, {:?}, {:?}) = {:?}, want {:?}",
                t.p1, t.p2, t.p3, result, t.want
            );
            // Test robust_sign(b,c,a) == robust_sign(a,b,c) for all a,b,c.
            let rotated = robust_sign(&t.p2, &t.p3, &t.p1);
            assert_eq!(
                rotated, result,
                "robust_sign({:?}, {:?}, {:?}) vs rotated robust_sign({:?}, {:?}, {:?}) = {:?}, want {:?}",
                t.p1, t.p2, t.p3, t.p2, t.p3, t.p1, rotated, result
            );
            // Test robust_sign(c,b,a) == -robust_sign(a,b,c) for all a,b,c.
            let want = match result {
                Direction::Clockwise => Direction::CounterClockwise,
                Direction::CounterClockwise => Direction::Clockwise,
                Direction::Indeterminate => Direction::Indeterminate,
            };
            let reversed = robust_sign(&t.p3, &t.p2, &t.p1);
            assert_eq!(
                reversed, want,
                "robust_sign({:?}, {:?}, {:?}) vs reversed robust_sign({:?}, {:?}, {:?}) = {:?}, want {:?}",
                t.p1, t.p2, t.p3, t.p3, t.p2, t.p1, reversed, want
            );
        }

        // None of the edge cases above should be indeterminate now that
        // exact_sign() is implemented.
        assert_ne!(
            robust_sign(&po_a(), &po_b(), &po_c()),
            Direction::Indeterminate
        );
        assert_ne!(
            robust_sign(&x1(), &x2(), &Point(x1().0 * -1.0)),
            Direction::Indeterminate
        );
        assert_ne!(
            robust_sign(&x3(), &x4(), &Point(x3().0 * -1.0)),
            Direction::Indeterminate
        );
        assert_ne!(
            robust_sign(&y1(), &y2(), &Point(y1().0 * -1.0)),
            Direction::Indeterminate
        );
    }

    #[test]
    fn test_stable_sign_underflow() {
        // Verify that stable_sign returns Indeterminate when its error
        // calculation underflows, while exact_sign (and therefore
        // robust_sign) still resolves the correct answer via exact
        // arithmetic.
        let a = pt(1., 1.9535722048627587e-90, 7.4882501322554515e-80);
        let b = pt(1., 9.6702373087191359e-127, 3.706704857169321e-116);
        let c = pt(1., 3.8163353663361477e-142, 1.4628419538608985e-131);

        assert_eq!(stable_sign(&a, &b, &c), Direction::Indeterminate);
        assert_eq!(exact_sign(&a, &b, &c, true), Direction::CounterClockwise);
        assert_eq!(robust_sign(&a, &b, &c), Direction::CounterClockwise);
    }

    // check_symbolic_sign verifies that expensive_sign resolves a, b, c
    // (which must be exactly coplanar with the origin, with a < b < c in
    // lexicographic order) to `expected` for every rotation of (a,b,c) and
    // to the opposite direction for every reversal. This is intended
    // specifically for checking the cases where symbolic perturbation is
    // needed to break ties.
    fn check_symbolic_sign(expected: Direction, a: Point, b: Point, c: Point) {
        assert!(a.0 < b.0);
        assert!(b.0 < c.0);
        assert_eq!(a.0.dot(&b.0.cross(&c.0)), 0.0);

        assert_eq!(expensive_sign(&a, &b, &c), expected);
        assert_eq!(expensive_sign(&b, &c, &a), expected);
        assert_eq!(expensive_sign(&c, &a, &b), expected);
        let reversed = flip(expected);
        assert_eq!(expensive_sign(&c, &b, &a), reversed);
        assert_eq!(expensive_sign(&b, &a, &c), reversed);
        assert_eq!(expensive_sign(&a, &c, &b), reversed);
    }

    /// Exercises every branch of symbolically_perturbed_sign(). Let M_1,
    /// M_2, ... be the sequence of submatrices whose determinant sign is
    /// tested by that function. Then the i-th case below is a 3x3 matrix M
    /// (with rows A, B, C) such that:
    ///
    ///    det(M) = 0
    ///    det(M_j) = 0 for j < i
    ///    det(M_i) != 0
    ///    A < B < C in lexicographic order.
    ///
    /// Reversing the sign of any "return" in symbolically_perturbed_sign()
    /// should make one of these cases fail.
    #[test]
    fn test_symbolic_perturbation_code_coverage() {
        // det(M_1) = b0*c1 - b1*c0
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(-3., -1., 0.),
            pt(-2., 1., 0.),
            pt(1., -2., 0.),
        );
        // det(M_2) = b2*c0 - b0*c2
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(-6., 3., 3.),
            pt(-4., 2., -1.),
            pt(-2., 1., 4.),
        );
        // det(M_3) = b1*c2 - b2*c1
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(0., -1., -1.),
            pt(0., 1., -2.),
            pt(0., 2., 1.),
        );
        // From this point onward, B or C must be zero, or B is proportional to C.
        // det(M_4) = c0*a1 - c1*a0
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(-1., 2., 7.),
            pt(2., 1., -4.),
            pt(4., 2., -8.),
        );
        // det(M_5) = c0
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(-4., -2., 7.),
            pt(2., 1., -4.),
            pt(4., 2., -8.),
        );
        // det(M_6) = -c1
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(0., -5., 7.),
            pt(0., -4., 8.),
            pt(0., -2., 4.),
        );
        // det(M_7) = c2*a0 - c0*a2
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(-5., -2., 7.),
            pt(0., 0., -2.),
            pt(0., 0., -1.),
        );
        // det(M_8) = c2
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(0., -2., 7.),
            pt(0., 0., 1.),
            pt(0., 0., 2.),
        );
        // From this point onward, C must be zero.
        // det(M_9) = a0*b1 - a1*b0
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(-3., 1., 7.),
            pt(-1., -4., 1.),
            pt(0., 0., 0.),
        );
        // det(M_10) = -b0
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(-6., -4., 7.),
            pt(-3., -2., 1.),
            pt(0., 0., 0.),
        );
        // det(M_11) = b1
        check_symbolic_sign(
            Direction::Clockwise,
            pt(0., -4., 7.),
            pt(0., -2., 1.),
            pt(0., 0., 0.),
        );
        // det(M_12) = a0
        check_symbolic_sign(
            Direction::Clockwise,
            pt(-1., -4., 5.),
            pt(0., 0., -3.),
            pt(0., 0., 0.),
        );
        // det(M_13) = 1
        check_symbolic_sign(
            Direction::CounterClockwise,
            pt(0., -4., 5.),
            pt(0., 0., -5.),
            pt(0., 0., 0.),
        );
    }
}
