//! This module defines [`EdgeCrosser`], which allows edges to be efficiently tested
//! for intersection with a given fixed edge `ab`. It is especially efficient when testing
//! for intersection with an edge chain connecting vertices `v₀`, `v₁`, `v₂`, ...

use super::{Crossing, SignedCrossing, robust_cross_prod, signed_vertex_crossing, vertex_crossing};
use crate::consts::*;
use crate::s2::point::Point;
use crate::s2::predicates::{Direction, robust_sign, triage_sign};

fn invert_dir(d: Direction) -> Direction {
    match d {
        Direction::CounterClockwise => Direction::Clockwise,
        Direction::Clockwise => Direction::CounterClockwise,
        Direction::Indeterminate => Direction::Indeterminate,
    }
}

/// `EdgeCrosser` allows edges to be efficiently tested for intersection with a
/// given fixed edge `ab`. It is especially efficient when testing for
/// intersection with an edge chain connecting vertices `v₀`, `v₁`, `v₂`, ...
///
/// # Examples
///
/// Testing against independent, disjoint edges:
/// ```rust
/// # use s2::point::Point;
/// # use s2::edge_crossings::{EdgeCrosser, Crossing};
/// fn count_intersections(a: &Point, b: &Point, edges: &[(Point, Point)]) -> usize {
///     let mut crosser = EdgeCrosser::new(a, b);
///     edges.iter()
///         .filter(|e| crosser.crossing_sign(&e.0, &e.1) >= Crossing::MaybeCross)
///         .count()
/// }
/// ```
///
/// Testing against a connected edge chain to avoid redundant vertex evaluations:
/// ```rust
/// # use s2::point::Point;
/// # use s2::edge_crossings::{EdgeCrosser, Crossing};
/// fn count_chain_intersections(a: &Point, b: &Point, chain: &[Point]) -> usize {
///     if chain.len() < 2 { return 0; }
///     
///     let mut count = 0;
///     let mut crosser = EdgeCrosser::new_with_c(a, b, &chain[0]);
///     for vertex in &chain[1..] {
///         if crosser.chain_crossing_sign(vertex) >= Crossing::MaybeCross {
///             count += 1;
///         }
///     }
///     count
/// }
/// ```
#[derive(Copy, Clone)]
pub struct EdgeCrosser {
    a: Point,
    b: Point,
    have_tangents: bool,
    a_tangent: Point,
    b_tangent: Point,
    c: Point,
    acb: Direction,
    bda: Direction,
}

impl EdgeCrosser {
    /// Creates an `EdgeCrosser` with the fixed edge `ab`.
    pub fn new(a: &Point, b: &Point) -> Self {
        debug_assert!(a.0.is_unit());
        debug_assert!(b.0.is_unit());
        EdgeCrosser {
            a: *a,
            b: *b,
            have_tangents: false,
            a_tangent: Point::default(),
            b_tangent: Point::default(),
            c: Point::default(),
            acb: Direction::Indeterminate,
            bda: Direction::Indeterminate,
        }
    }

    /// Creates an `EdgeCrosser` with the fixed edge `ab` and the first vertex `c`.
    pub fn new_with_c(a: &Point, b: &Point, c: &Point) -> Self {
        let mut crosser = Self::new(a, b);
        crosser.restart_at(c);
        crosser
    }

    /// Sets the current vertex of the chain to `c`.
    pub fn restart_at(&mut self, c: &Point) {
        debug_assert!(c.0.is_unit());
        self.c = *c;
        self.acb = invert_dir(triage_sign(&self.a, &self.b, &self.c));
    }

    /// Reports whether the edge `ab` intersects the edge `cd`.
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
    /// that even if you have 4 points `A, B, C, D` that lie exactly in a line
    /// (say, around the equator), `C` and `D` will be treated as being slightly to
    /// one side or the other of `AB`. This is done in a way such that the
    /// results are always consistent with [`robust_sign`].
    ///
    /// Note that if you want to check an edge against a chain of other edges,
    /// it is slightly more efficient to use the [`EdgeCrosser::chain_crossing_sign`] method.
    pub fn crossing_sign(&mut self, c: &Point, d: &Point) -> Crossing {
        debug_assert!(c.0.is_unit());
        debug_assert!(d.0.is_unit());
        if c != &self.c {
            self.restart_at(c);
        }
        self.chain_crossing_sign(d)
    }

    // Note: In C++ S2Geometry, these single-argument chain methods are function overloads
    // (e.g. CrossingSign(c, d) and CrossingSign(d)). Since Rust does not support method
    // overloading, we prefix the chain-based equivalents with `chain_`.

    /// Behaves like [`EdgeCrosser::crossing_sign`] but uses the last vertex
    /// passed to one of the crossing methods (or [`EdgeCrosser::restart_at`]) as the first vertex of the current edge.
    pub fn chain_crossing_sign(&mut self, d: &Point) -> Crossing {
        debug_assert!(d.0.is_unit());
        let bda = triage_sign(&self.a, &self.b, d);
        if self.acb == invert_dir(bda) && bda != Direction::Indeterminate {
            self.c = *d;
            self.acb = invert_dir(bda);
            return Crossing::DoNotCross;
        }
        self.bda = bda;
        self.crossing_sign_internal(d)
    }

    /// This method extends the concept of a "crossing" to the case where `ab`
    /// and `cd` have a vertex in common. The two edges may or may not cross,
    /// according to the rules defined in [`vertex_crossing`]. The rules
    /// are designed so that point containment tests can be implemented simply
    /// by counting edge crossings. Similarly, determining whether one edge
    /// chain crosses another edge chain can be implemented by counting.
    ///
    /// Returns true if `crossing_sign(c, d)` is [`Crossing::Cross`], or `ab` and `cd` share a vertex
    /// and `vertex_crossing(a, b, c, d)` returns true.
    pub fn edge_or_vertex_crossing(&mut self, c: &Point, d: &Point) -> bool {
        debug_assert!(c.0.is_unit());
        debug_assert!(d.0.is_unit());
        if c != &self.c {
            self.restart_at(c);
        }
        self.chain_edge_or_vertex_crossing(d)
    }

    pub fn chain_edge_or_vertex_crossing(&mut self, d: &Point) -> bool {
        debug_assert!(d.0.is_unit());
        let c = self.c; // copy since chain_crossing_sign clobbers self.c
        let crossing = self.chain_crossing_sign(d);
        match crossing {
            Crossing::DoNotCross => false,
            Crossing::Cross => true,
            Crossing::MaybeCross => vertex_crossing(&self.a, &self.b, &c, d),
        }
    }

    /// Like [`EdgeCrosser::edge_or_vertex_crossing`], but returns [`SignedCrossing::LeftToRight`] if `ab` crosses `cd` from left to
    /// right, [`SignedCrossing::RightToLeft`] if `ab` crosses `cd` from right to left, and [`SignedCrossing::None`] otherwise. This
    /// implies that if `cd` bounds some region according to the "interior is on
    /// the left" rule, this function returns [`SignedCrossing::LeftToRight`] when `ab` exits the region and [`SignedCrossing::RightToLeft`]
    /// when `ab` enters.
    ///
    /// This method allows computing the change in winding number from point `a` to
    /// point `b` by summing the signed edge crossings of `ab` with the edges of the
    /// loop(s) used to define the winding number.
    pub fn signed_edge_or_vertex_crossing(&mut self, c: &Point, d: &Point) -> SignedCrossing {
        debug_assert!(c.0.is_unit());
        debug_assert!(d.0.is_unit());
        if c != &self.c {
            self.restart_at(c);
        }
        self.chain_signed_edge_or_vertex_crossing(d)
    }

    pub fn chain_signed_edge_or_vertex_crossing(&mut self, d: &Point) -> SignedCrossing {
        debug_assert!(d.0.is_unit());
        let c = self.c;
        let crossing = self.chain_crossing_sign(d);
        match crossing {
            Crossing::DoNotCross => SignedCrossing::None,
            Crossing::Cross => match self.last_interior_crossing_sign() {
                Direction::CounterClockwise => SignedCrossing::RightToLeft,
                Direction::Clockwise => SignedCrossing::LeftToRight,
                Direction::Indeterminate => SignedCrossing::None,
            },
            Crossing::MaybeCross => signed_vertex_crossing(&self.a, &self.b, &c, d),
        }
    }

    /// If the preceding call to [`EdgeCrosser::crossing_sign`] returned [`Crossing::Cross`], this method returns the direction
    /// of the crossing.
    pub fn last_interior_crossing_sign(&self) -> Direction {
        self.acb
    }

    fn crossing_sign_internal(&mut self, d: &Point) -> Crossing {
        let result = self.crossing_sign_internal2(d);
        self.c = *d;
        self.acb = invert_dir(self.bda);
        result
    }

    fn crossing_sign_internal2(&mut self, d: &Point) -> Crossing {
        if !self.have_tangents {
            let norm = robust_cross_prod(&self.a, &self.b);
            self.a_tangent = Point(self.a.0.cross(&norm.0));
            self.b_tangent = Point(norm.0.cross(&self.b.0));
            self.have_tangents = true;
        }

        let error_bound = (1.5 + 1.0 / 3.0f64.sqrt()) * DBL_EPSILON;
        if (self.c.0.dot(&self.a_tangent.0) > error_bound
            && d.0.dot(&self.a_tangent.0) > error_bound)
            || (self.c.0.dot(&self.b_tangent.0) > error_bound
                && d.0.dot(&self.b_tangent.0) > error_bound)
        {
            return Crossing::DoNotCross;
        }

        if self.a == self.c || self.a == *d || self.b == self.c || self.b == *d {
            return Crossing::MaybeCross;
        }

        if self.a == self.b || self.c == *d {
            return Crossing::DoNotCross;
        }

        if self.acb == Direction::Indeterminate {
            self.acb = invert_dir(robust_sign(&self.a, &self.b, &self.c));
        }
        if self.bda == Direction::Indeterminate {
            self.bda = robust_sign(&self.a, &self.b, d);
        }
        if self.bda != self.acb {
            return Crossing::DoNotCross;
        }

        let cbd = invert_dir(robust_sign(&self.c, d, &self.b));
        if cbd != self.acb {
            return Crossing::DoNotCross;
        }

        let dac = robust_sign(&self.c, d, &self.a);
        if dac != self.acb {
            Crossing::DoNotCross
        } else {
            Crossing::Cross
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r3::vector::Vector;
    use crate::s2::point::Point;

    #[test]
    #[should_panic]
    fn test_invalid_default_points() {
        let p = Point::default();
        let _crosser = EdgeCrosser::new(&p, &p);
    }

    #[test]
    #[should_panic]
    fn test_invalid_nan_points() {
        let nan = std::f64::NAN;
        let p = Point(Vector::new(nan, nan, nan));
        let _crosser = EdgeCrosser::new(&p, &p);
    }

    fn test_crossing(
        a: &Point,
        b: &Point,
        c: &Point,
        d: &Point,
        crossing_sign: i32,
        signed_crossing_sign: i32,
    ) {
        let mut crossing_sign = crossing_sign;
        if a == c || a == d || b == c || b == d {
            crossing_sign = 0;
        }

        let expected_cross = match crossing_sign {
            1 => Crossing::Cross,
            0 => Crossing::MaybeCross,
            -1 => Crossing::DoNotCross,
            _ => unreachable!(),
        };

        let edge_or_vertex = signed_crossing_sign != 0;

        let expected_signed = match signed_crossing_sign {
            1 => SignedCrossing::RightToLeft,
            0 => SignedCrossing::None,
            -1 => SignedCrossing::LeftToRight,
            _ => unreachable!(),
        };

        assert_eq!(
            crate::s2::edge_crossings::crossing_sign(a, b, c, d),
            expected_cross
        );
        if a == c || a == d || b == c || b == d {
            assert_eq!(
                crate::s2::edge_crossings::edge_or_vertex_crossing(a, b, c, d),
                edge_or_vertex
            );
            assert_eq!(
                crate::s2::edge_crossings::signed_vertex_crossing(a, b, c, d) as i32,
                expected_signed as i32
            );
        }

        let mut crosser = EdgeCrosser::new_with_c(a, b, c);
        assert_eq!(crosser.chain_crossing_sign(d), expected_cross);
        assert_eq!(crosser.chain_crossing_sign(c), expected_cross);

        if a == c || a == d || b == c || b == d {
            crosser.restart_at(c);
            assert_eq!(crosser.chain_edge_or_vertex_crossing(d), edge_or_vertex);
            assert_eq!(crosser.chain_edge_or_vertex_crossing(c), edge_or_vertex);

            crosser.restart_at(c);
            assert_eq!(
                crosser.chain_signed_edge_or_vertex_crossing(d),
                expected_signed
            );
            assert_eq!(
                crosser
                    .chain_signed_edge_or_vertex_crossing(c)
                    .winding_delta(),
                -expected_signed.winding_delta()
            );
        }
    }

    fn test_crossings_permutations(
        a: Point,
        b: Point,
        c: Point,
        d: Point,
        crossing_sign: i32,
        signed_crossing_sign: i32,
    ) {
        let a = Point(a.0.normalize());
        let b = Point(b.0.normalize());
        let c = Point(c.0.normalize());
        let d = Point(d.0.normalize());
        test_crossing(&a, &b, &c, &d, crossing_sign, signed_crossing_sign);
        test_crossing(&b, &a, &c, &d, crossing_sign, -signed_crossing_sign);
        test_crossing(&a, &b, &d, &c, crossing_sign, -signed_crossing_sign);
        test_crossing(&b, &a, &d, &c, crossing_sign, signed_crossing_sign);
        test_crossing(&a, &a, &c, &d, -1, 0);
        test_crossing(&a, &b, &c, &c, -1, 0);
        test_crossing(&a, &a, &c, &c, -1, 0);
        test_crossing(&a, &b, &a, &b, 0, 1);
        if crossing_sign == 0 {
            test_crossing(&c, &d, &a, &b, crossing_sign, 0);
        } else {
            test_crossing(&c, &d, &a, &b, crossing_sign, -signed_crossing_sign);
        }
    }

    #[test]
    fn test_crossings() {
        test_crossings_permutations(
            Point(Vector::new(1.0, 2.0, 1.0)),
            Point(Vector::new(1.0, -3.0, 0.5)),
            Point(Vector::new(1.0, -0.5, -3.0)),
            Point(Vector::new(0.1, 0.5, 3.0)),
            1,
            1,
        );

        test_crossings_permutations(
            Point(Vector::new(1.0, 2.0, 1.0)),
            Point(Vector::new(1.0, -3.0, 0.5)),
            Point(Vector::new(-1.0, 0.5, 3.0)),
            Point(Vector::new(-0.1, -0.5, -3.0)),
            -1,
            0,
        );

        test_crossings_permutations(
            Point(Vector::new(0.0, 0.0, -1.0)),
            Point(Vector::new(0.0, 1.0, 0.0)),
            Point(Vector::new(0.0, 0.0, 1.0)),
            Point(Vector::new(0.0, 1.0, 1.0)),
            -1,
            0,
        );

        test_crossings_permutations(
            Point(Vector::new(1.0, 0.0, 0.0)),
            Point(Vector::new(0.0, 1.0, 0.0)),
            Point(Vector::new(1.0, -0.1, 1.0)),
            Point(Vector::new(1.0, 1.0, -0.1)),
            1,
            1,
        );

        test_crossings_permutations(
            Point(Vector::new(7.0, -2.0, 3.0)),
            Point(Vector::new(2.0, 3.0, 4.0)),
            Point(Vector::new(2.0, 3.0, 4.0)),
            Point(Vector::new(-1.0, 2.0, 5.0)),
            0,
            -1,
        );

        test_crossings_permutations(
            Point(Vector::new(1.0, 1.0, 1.0)),
            Point(Vector::new(1.0, 1.0 + std::f64::EPSILON, -1.0)),
            Point(Vector::new(1.0, -1.0, 0.0)),
            Point(Vector::new(1.0, 1.0, 0.0)),
            -1,
            0,
        );
    }

    #[test]
    fn test_collinear_edges_that_dont_touch() {
        let a = Point(Vector::new(1.0, 1.0, 1.0).normalize());
        let b = Point(Vector::new(1.0, -1.0, 1.0).normalize());
        let c = Point(Vector::new(-1.0, -1.0, 1.0).normalize());
        let d = Point(Vector::new(-1.0, 1.0, 1.0).normalize());
        let mut crosser = EdgeCrosser::new_with_c(&a, &b, &c);
        assert_eq!(crosser.chain_crossing_sign(&d), Crossing::DoNotCross);
    }

    #[test]
    fn test_coincident_zero_length_edges_that_dont_touch() {
        let a = Point(Vector::new(1.0, 1.0, 1.0).normalize());
        let mut crosser = EdgeCrosser::new_with_c(&a, &a, &a);
        assert_eq!(crosser.chain_crossing_sign(&a), Crossing::MaybeCross);
    }
}
