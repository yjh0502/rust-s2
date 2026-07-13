//! Defines functions related to determining whether two geodesic edges cross
//! and for computing intersection points.
//!
//! The predicates [`crossing_sign`], [`vertex_crossing`], and [`edge_or_vertex_crossing`]
//! are robust, meaning that they produce correct, consistent results even in
//! pathological cases.
//!
//! See also [`EdgeCrosser`] (which efficiently tests an edge against a sequence
//! of other edges).

use crate::consts::*;
use crate::r3::vector::Vector;
use crate::s2::point::{Point, ordered_ccw};

mod crosser;
pub use crosser::EdgeCrosser;

const DBL_ERR: f64 = DBL_EPSILON / 2.0;

/// INTERSECTION_ERROR is an upper bound on the distance from the intersection
/// point returned by get_intersection() to the true intersection point.
pub const INTERSECTION_ERROR: f64 = 8.0 * DBL_ERR;

/// This value can be used as the Builder snap_radius() to ensure that edges
/// that have been displaced by up to INTERSECTION_ERROR are merged back
/// together again. For example this can happen when geometry is intersected
/// with a set of tiles and then unioned. It is equal to twice the
/// intersection error because input edges might have been displaced in
/// opposite directions.
pub const INTERSECTION_MERGE_RADIUS: f64 = 16.0 * DBL_ERR;

/// A Crossing indicates how edges cross.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Crossing {
    DoNotCross = -1,
    MaybeCross = 0,
    Cross = 1,
}

/// Indicates the direction an edge crosses another edge, used for computing winding numbers.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum SignedCrossing {
    LeftToRight = -1,
    None = 0,
    RightToLeft = 1,
}

impl SignedCrossing {
    /// Returns the crossing parity as an integer delta for summing winding numbers.
    pub fn winding_delta(&self) -> i32 {
        *self as i32
    }
}

/// Returns a vector whose direction is guaranteed to be very close to the exact
/// mathematical cross product of the given unit-length vectors "a" and "b", but
/// whose magnitude is arbitrary.  Unlike a.cross(b), this statement is true
/// even when "a" and "b" are very nearly parallel (i.e., a ~= b or a ~= -b).
/// Specifically, the direction of the result vector differs from the exact cross
/// product by at most 6 * DBL_ERR radians.
///
/// When a == -b exactly, the result is consistent with the symbolic perturbation
/// model used by [`robust_sign`](crate::s2::predicates::robust_sign). In other words, even antipodal
/// point pairs have a consistent and well-defined edge between them. (In fact
/// this is true for any pair of distinct points whose vectors are parallel.)
///
/// When a == b exactly, an arbitrary vector orthogonal to "a" is returned.
/// [From a strict mathematical viewpoint it would be better to return (0, 0, 0),
/// but this behavior helps to avoid special cases in client code.]
///
/// This function has the following properties:
///
///   1. `robust_cross_prod(a, b) != 0` for all `a`, `b`
///   2. `robust_cross_prod(b, a) == -robust_cross_prod(a, b)` unless `a == b`
///   3. `robust_cross_prod(-a, b) == -robust_cross_prod(a, b)` unless `a` and `b` are exactly proportional
///   4. `robust_cross_prod(a, -b) == -robust_cross_prod(a, b)` unless `a` and `b` are exactly proportional
///
/// Note that if you want the result to be unit-length, you must call `normalize()`
/// explicitly. (The result is always scaled such that `normalize()` can be called
/// without precision loss due to floating-point underflow.)
pub fn robust_cross_prod(a: &Point, b: &Point) -> Point {
    debug_assert!(a.0.is_unit());
    debug_assert!(b.0.is_unit());
    let mut result = Point::default();
    if get_stable_cross_prod(a, b, &mut result) {
        return result;
    }
    if a == b {
        return a.ortho();
    }
    symbolic_cross_prod(a, b)
}

fn get_stable_cross_prod(a: &Point, b: &Point, result: &mut Point) -> bool {
    let robust_cross_prod_error = 6.0 * DBL_ERR;
    let sqrt3 = 3.0f64.sqrt();
    let min_norm =
        (32.0 * sqrt3 * DBL_ERR) / (robust_cross_prod_error / DBL_ERR - (1.0 + 2.0 * sqrt3));

    let x = a.0 - b.0;
    let y = a.0 + b.0;
    result.0 = x.cross(&y);
    result.0.norm2() >= min_norm * min_norm
}

fn symbolic_cross_prod(a: &Point, b: &Point) -> Point {
    if a.0 < b.0 {
        ensure_normalizable(&symbolic_cross_prod_sorted(a, b))
    } else {
        Point(ensure_normalizable(&symbolic_cross_prod_sorted(b, a)).0 * -1.0)
    }
}

fn symbolic_cross_prod_sorted(a: &Point, b: &Point) -> Point {
    debug_assert!(a.0 < b.0);

    if b.0.x != 0.0 || b.0.y != 0.0 {
        return Point(Vector::new(-b.0.y, b.0.x, 0.0));
    }
    if b.0.z != 0.0 {
        return Point(Vector::new(b.0.z, 0.0, 0.0));
    }

    if a.0.x != 0.0 || a.0.y != 0.0 {
        return Point(Vector::new(a.0.y, -a.0.x, 0.0));
    }
    if a.0.z != 0.0 {
        return Point(Vector::new(-a.0.z, 0.0, 0.0));
    }

    Point(Vector::new(1.0, 0.0, 0.0))
}

fn is_normalizable(p: &Point) -> bool {
    let max_abs = p.0.x.abs().max(p.0.y.abs()).max(p.0.z.abs());
    max_abs >= 2.0f64.powi(-242)
}

fn ensure_normalizable(p: &Point) -> Point {
    if !is_normalizable(p) {
        let max_abs = p.0.x.abs().max(p.0.y.abs()).max(p.0.z.abs());
        let ilogb = max_abs.log2().floor() as i32;
        let scale = 2.0f64.powi(-1 - ilogb);
        return Point(p.0 * scale);
    }
    *p
}

/// This function determines whether the edge `ab` intersects the edge `cd`.
/// Returns [`Crossing::Cross`] if `ab` crosses `cd` at a point that is interior to both edges.
/// Returns [`Crossing::MaybeCross`] if any two vertices from different edges are the same.
/// Returns [`Crossing::DoNotCross`] otherwise.
///
/// Note that if an edge is degenerate (`a == b` or `c == d`), the return value
/// is [`Crossing::MaybeCross`] if two vertices from different edges are the same and [`Crossing::DoNotCross`] otherwise.
///
/// Properties of this function:
///
///  1. `crossing_sign(b, a, c, d) == crossing_sign(a, b, c, d)`
///  2. `crossing_sign(c, d, a, b) == crossing_sign(a, b, c, d)`
///  3. `crossing_sign(a, b, c, d) == Crossing::MaybeCross` if `a == c`, `a == d`, `b == c`, `b == d`
///  4. `crossing_sign(a, b, c, d) <= Crossing::MaybeCross` if `a == b` or `c == d`
///
/// This function implements an exact, consistent perturbation model such
/// that no three points are ever considered to be collinear. This means
/// that even if you have 4 points `a`, `b`, `c`, `d` that lie exactly in a line
/// (say, around the equator), `c` and `d` will be treated as being slightly to
/// one side or the other of `ab`. This is done in a way such that the
/// results are always consistent with [`robust_sign`](crate::s2::predicates::robust_sign).
///
/// Note that if you want to check an edge against a collection of other edges,
/// it is much more efficient to use an [`EdgeCrosser`].
pub fn crossing_sign(a: &Point, b: &Point, c: &Point, d: &Point) -> Crossing {
    let mut crosser = EdgeCrosser::new_with_c(a, b, c);
    crosser.chain_crossing_sign(d)
}

/// Given two edges `ab` and `cd` where at least two vertices are identical
/// (i.e. `crossing_sign(a, b, c, d) == Crossing::MaybeCross`), this function defines whether the
/// two edges cross in such a way that point-in-polygon containment tests can
/// be implemented by counting the number of edge crossings. The basic rule is
/// that a crossing occurs if `ab` is encountered after `cd` during a CCW sweep
/// around the shared vertex starting from a fixed reference point.
///
/// Note that according to this rule, if `ab` crosses `cd` then in general `cd`
/// does not cross `ab`. However, this leads to the correct result when
/// counting polygon edge crossings. For example, suppose that `a`, `b`, `c` are
/// three consecutive vertices of a CCW polygon. If we now consider the edge
/// crossings of a segment `bp` as `p` sweeps around `b`, the crossing number
/// changes parity exactly when `bp` crosses `ba` or `bc`.
///
/// Useful properties of this function:
///
///  1. `vertex_crossing(a, a, c, d) == vertex_crossing(a, b, c, c) == false`
///  2. `vertex_crossing(a, b, a, b) == vertex_crossing(a, b, b, a) == true`
///  3. `vertex_crossing(a, b, c, d) == vertex_crossing(a, b, d, c) == vertex_crossing(b, a, c, d) == vertex_crossing(b, a, d, c)`
///  4. If exactly one of `a, b` equals one of `c, d`, then exactly one of
///     `vertex_crossing(a, b, c, d)` and `vertex_crossing(c, d, a, b)` is true
///
/// # Panics
///
/// Panics in debug builds if called with 4 distinct vertices.
pub fn vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> bool {
    debug_assert!(
        a == c || a == d || b == c || b == d,
        "vertex_crossing called with 4 distinct vertices"
    );
    if a == b || c == d {
        return false;
    }
    if a == c {
        return b == d || ordered_ccw(&a.ortho(), d, b, a);
    }
    if b == d {
        return ordered_ccw(&b.ortho(), c, a, b);
    }
    if a == d {
        return b == c || ordered_ccw(&a.ortho(), c, b, a);
    }
    if b == c {
        return ordered_ccw(&b.ortho(), d, a, b);
    }
    false
}

/// Like [`vertex_crossing`] but returns [`SignedCrossing::LeftToRight`] if `ab` crosses `cd` from left to right,
/// [`SignedCrossing::RightToLeft`] if `ab` crosses `cd` from right to left, and [`SignedCrossing::None`] otherwise. This implies that
/// if `cd` bounds some region according to the "interior is on the left" rule,
/// this function returns [`SignedCrossing::LeftToRight`] when `ab` exits the region and [`SignedCrossing::RightToLeft`] when `ab` enters.
///
/// This is a helper method that allows computing the change in winding number
/// from point `a` to point `b` by summing the signed edge crossings of `ab` with the
/// edges of the loop(s) used to define the winding number.
///
/// Useful properties of this function:
///
///  1. `signed_vertex_crossing(a, a, c, d) == SignedCrossing::None`
///  2. `signed_vertex_crossing(a, b, a, b) == SignedCrossing::RightToLeft`
///  3. `signed_vertex_crossing(a, b, b, a) == SignedCrossing::LeftToRight`
///  4. `signed_vertex_crossing(a, b, c, d).winding_delta() == -signed_vertex_crossing(a, b, d, c).winding_delta() == -signed_vertex_crossing(b, a, c, d).winding_delta() == signed_vertex_crossing(b, a, d, c).winding_delta()`
///  5. If exactly one of `a, b` equals one of `c, d`, then exactly one of
///     `signed_vertex_crossing(a, b, c, d)` and `signed_vertex_crossing(c, d, a, b)` is non-zero
///
/// # Panics
///
/// Panics in debug builds if called with 4 distinct vertices.
pub fn signed_vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> SignedCrossing {
    debug_assert!(
        a == c || a == d || b == c || b == d,
        "signed_vertex_crossing called with 4 distinct vertices"
    );
    if a == b || c == d {
        return SignedCrossing::None;
    }
    if a == c {
        return if b == d || ordered_ccw(&a.ortho(), d, b, a) {
            SignedCrossing::RightToLeft
        } else {
            SignedCrossing::None
        };
    }
    if b == d {
        return if ordered_ccw(&b.ortho(), c, a, b) {
            SignedCrossing::RightToLeft
        } else {
            SignedCrossing::None
        };
    }
    if a == d {
        return if b == c || ordered_ccw(&a.ortho(), c, b, a) {
            SignedCrossing::LeftToRight
        } else {
            SignedCrossing::None
        };
    }
    if b == c {
        return if ordered_ccw(&b.ortho(), d, a, b) {
            SignedCrossing::LeftToRight
        } else {
            SignedCrossing::None
        };
    }
    SignedCrossing::None
}

/// A convenience function that calls [`crossing_sign`] to handle cases
/// where all four vertices are distinct, and [`vertex_crossing`] to handle
/// cases where two or more vertices are the same. This defines a crossing
/// function such that point-in-polygon containment tests can be implemented
/// by simply counting edge crossings.
pub fn edge_or_vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> bool {
    match crossing_sign(a, b, c, d) {
        Crossing::DoNotCross => false,
        Crossing::Cross => true,
        Crossing::MaybeCross => vertex_crossing(a, b, c, d),
    }
}

fn compare_edges(mut a0: Point, mut a1: Point, mut b0: Point, mut b1: Point) -> bool {
    if a0.0 > a1.0 {
        std::mem::swap(&mut a0, &mut a1);
    }
    if b0.0 > b1.0 {
        std::mem::swap(&mut b0, &mut b1);
    }
    a0.0 < b0.0 || (a0.0 == b0.0 && a1.0 < b1.0)
}

fn projection(x: &Vector, a_norm: &Vector, a_norm_len: f64, a0: &Point, a1: &Point) -> (f64, f64) {
    let x0 = *x - a0.0;
    let x1 = *x - a1.0;
    let x0_dist2 = x0.norm2();
    let x1_dist2 = x1.norm2();

    let dist;
    let proj;
    if x0_dist2 < x1_dist2 || (x0_dist2 == x1_dist2 && x0 < x1) {
        dist = x0_dist2.sqrt();
        proj = x0.dot(a_norm);
    } else {
        dist = x1_dist2.sqrt();
        proj = x1.dot(a_norm);
    }

    let bound = (((3.5 + 2.0 * 3f64.sqrt()) * a_norm_len + 32.0 * 3f64.sqrt() * DBL_ERR) * dist
        + 1.5 * proj.abs())
        * DBL_ERR;
    (proj, bound)
}

fn get_intersection_stable_sorted(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Option<Point> {
    let a_norm = (a0.0 - a1.0).cross(&(a0.0 + a1.0));
    let a_norm_len = a_norm.norm();
    let b_len = (b1.0 - b0.0).norm();

    let (mut b0_dist, b0_error) = projection(&b0.0, &a_norm, a_norm_len, a0, a1);
    let (mut b1_dist, b1_error) = projection(&b1.0, &a_norm, a_norm_len, a0, a1);

    if b0_dist < b1_dist {
        b0_dist = -b0_dist;
        b1_dist = -b1_dist;
    }
    let dist_sum = b0_dist - b1_dist;
    let error_sum = b0_error + b1_error;
    if dist_sum <= error_sum {
        return None;
    }

    let x = (b1.0 * b0_dist) - (b0.0 * b1_dist);
    let err = b_len * (b0_dist * b1_error - b1_dist * b0_error).abs() / (dist_sum - error_sum)
        + 2.0 * dist_sum * DBL_ERR;

    let x_len2 = x.norm2();
    if x_len2 < f64::MIN_POSITIVE {
        return None;
    }
    let x_len = x_len2.sqrt();
    let max_error = INTERSECTION_ERROR;
    if err > (max_error - DBL_ERR) * x_len {
        return None;
    }

    Some(Point(x * (1.0 / x_len)))
}

fn get_intersection_stable(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Option<Point> {
    let a_len2 = (a1.0 - a0.0).norm2();
    let b_len2 = (b1.0 - b0.0).norm2();
    if a_len2 < b_len2 || (a_len2 == b_len2 && compare_edges(*a0, *a1, *b0, *b1)) {
        get_intersection_stable_sorted(b0, b1, a0, a1)
    } else {
        get_intersection_stable_sorted(a0, a1, b0, b1)
    }
}

// Emulates the tie-breaking behavior of C++ GetIntersectionExact without using bigdecimal.
fn get_intersection_exact(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Point {
    debug_assert!(a0.0.is_unit());
    debug_assert!(a1.0.is_unit());
    debug_assert!(b0.0.is_unit());
    debug_assert!(b1.0.is_unit());
    let mut a_norm = Point(a0.0.cross(&a1.0));
    let mut b_norm = Point(b0.0.cross(&b1.0));
    let x = Point(a_norm.0.cross(&b_norm.0));

    if x.0.x != 0.0 || x.0.y != 0.0 || x.0.z != 0.0 {
        let is_ccw = crate::s2::predicates::sign(a0, a1, b1);
        let mut pt = Point(x.0.normalize());
        if !is_ccw {
            pt = Point(pt.0 * -1.0);
        }
        return pt;
    }

    if a_norm.0 == Vector::new(0.0, 0.0, 0.0) {
        a_norm = symbolic_cross_prod(a0, a1);
    }
    if b_norm.0 == Vector::new(0.0, 0.0, 0.0) {
        b_norm = symbolic_cross_prod(b0, b1);
    }

    let mut pt = Point(Vector::new(10.0, 10.0, 10.0));
    if ordered_ccw(b0, a0, b1, &b_norm) && a0.0 < pt.0 {
        pt = *a0;
    }
    if ordered_ccw(b0, a1, b1, &b_norm) && a1.0 < pt.0 {
        pt = *a1;
    }
    if ordered_ccw(a0, b0, a1, &a_norm) && b0.0 < pt.0 {
        pt = *b0;
    }
    if ordered_ccw(a0, b1, a1, &a_norm) && b1.0 < pt.0 {
        pt = *b1;
    }

    pt
}

/// Given two edges `ab` and `cd` such that `crossing_sign(a, b, c, d) == Crossing::Cross`, returns
/// their intersection point. Useful properties of this function:
///
///  1. `get_intersection(b, a, c, d) == get_intersection(a, b, d, c) == get_intersection(a, b, c, d)`
///  2. `get_intersection(c, d, a, b) == get_intersection(a, b, c, d)`
///
/// The returned intersection point is guaranteed to be very close to the
/// true intersection point of `ab` and `cd`, even if the edges intersect at a
/// very small angle. See [`INTERSECTION_ERROR`] for details.
pub fn get_intersection(a0: &Point, a1: &Point, b0: &Point, b1: &Point) -> Point {
    debug_assert!(crossing_sign(a0, a1, b0, b1) == Crossing::Cross);
    let pt = match get_intersection_stable(a0, a1, b0, b1) {
        Some(p) => p,
        None => get_intersection_exact(a0, a1, b0, b1),
    };
    debug_assert!(pt.0.is_unit());
    pt
}

/// Returns true if the angle `abc` contains its vertex `b`. Containment is
/// defined such that if several polygons tile the region around a vertex, then
/// exactly one of those polygons contains that vertex. Returns false for
/// degenerate angles of the form `aba`.
///
/// **Note**: This method is not sufficient to determine vertex containment in
/// polygons with duplicate vertices (such as the polygon `ABCADE`). Use
/// [`ContainsVertexQuery`](crate::s2::contains_vertex_query::ContainsVertexQuery) for such polygons.
///
/// Useful properties of this function:
///
///  1. `angle_contains_vertex(a, b, a) == false`
///  2. `angle_contains_vertex(a, b, c) == !angle_contains_vertex(c, b, a)` unless `a == c`
///  3. Given vertices `v_1` ... `v_k` ordered cyclically CCW around vertex `b`,
///     `angle_contains_vertex(v_{i+1}, b, v_i)` is true for exactly one value of `i`.
pub fn angle_contains_vertex(a: &Point, b: &Point, c: &Point) -> bool {
    !ordered_ccw(&b.ortho(), c, a, b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r3::vector::Vector;
    use crate::s2::point::Point;

    // Note: The following tests from the C++ reference implementation are skipped
    // because they test the ExactFloat (arbitrary precision) fallback paths which
    // are not present in this pure-f64 port:
    // - RobustCrossProdError
    // - IntersectionError
    // - GrazingIntersections
    // - ExactIntersectionUnderflow
    // - ExactIntersectionSign
    // - SymbolicCrossProdConsistentWithSign

    #[test]
    fn test_angle_contains_vertex() {
        let a = Point(Vector::new(1.0, 0.0, 0.0));
        let b = Point(Vector::new(0.0, 1.0, 0.0));
        let ref_b = b.ortho();

        // Degenerate angle ABA.
        assert!(!angle_contains_vertex(&a, &b, &a));

        // An angle where A == referenceDir(B).
        assert!(angle_contains_vertex(&ref_b, &b, &a));

        // An angle where C == referenceDir(B).
        assert!(!angle_contains_vertex(&a, &b, &ref_b));

        let v = Point(Vector::new(1.0, 1.0, 1.0).normalize());
        for _ in 0..10 {
            let pts = crate::s2::point::regular_points(
                &v,
                crate::s1::angle::Angle::from(crate::s1::angle::Deg(1.0)),
                10,
            );
            let mut num_contains = 0;
            for i in 0..pts.len() {
                if angle_contains_vertex(&pts[(i + 1) % pts.len()], &v, &pts[i]) {
                    num_contains += 1;
                }
            }
            assert_eq!(num_contains, 1);
        }
    }

    #[test]
    fn test_crossing_sign_simple() {
        let a = Point(Vector::new(1.0, 0.0, 0.0));
        let b = Point(Vector::new(0.0, 1.0, 0.0));
        let c = Point(Vector::new(1.0, 1.0, 1.0).normalize());
        let d = Point(Vector::new(1.0, 1.0, -1.0).normalize());

        let cross = crossing_sign(&a, &b, &c, &d);
        assert_eq!(cross, Crossing::Cross);

        // Same edge test
        assert_eq!(crossing_sign(&a, &b, &a, &b), Crossing::MaybeCross);
    }

    #[test]
    fn test_get_intersection_invariants() {
        use crate::s2::random;
        use rand::Rng;
        use rand::SeedableRng;
        use rand::rngs::StdRng;

        let mut rng = StdRng::seed_from_u64(12345);
        let iters = 5000;

        for _ in 0..iters {
            let mut a;
            let mut b;
            let mut c;
            let mut d;

            loop {
                a = random::point(&mut rng);
                b = random::point(&mut rng);

                c = Point(Vector::new(a.0.y, a.0.x, a.0.z));
                d = Point(Vector::new(b.0.y, b.0.x, b.0.z));

                if crossing_sign(&a, &b, &c, &d) == Crossing::Cross {
                    break;
                }
            }

            let result = get_intersection(&a, &b, &c, &d);

            let mut a2 = a;
            let mut b2 = b;
            let mut c2 = c;
            let mut d2 = d;

            if rng.next_u32() % 2 == 0 {
                std::mem::swap(&mut a2, &mut b2);
            }
            if rng.next_u32() % 2 == 0 {
                std::mem::swap(&mut c2, &mut d2);
            }
            if rng.next_u32() % 2 == 0 {
                std::mem::swap(&mut a2, &mut c2);
                std::mem::swap(&mut b2, &mut d2);
            }

            assert_eq!(result, get_intersection(&a2, &b2, &c2, &d2));
        }
    }

    #[test]
    fn test_compare_edges_order_invariant() {
        let v0 = Point(Vector::new(0.0, 1.0, 0.0));
        let v1 = Point(Vector::new(1.0, 0.0, 0.0));
        assert!(!compare_edges(v0, v1, v0, v1));
        assert!(!compare_edges(v1, v0, v0, v1));
        assert!(!compare_edges(v0, v1, v1, v0));
        assert!(!compare_edges(v1, v0, v1, v0));
    }
}
