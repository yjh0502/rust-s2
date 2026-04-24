/*
Copyright 2026 The rust-s2 Authors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
*/

//! Robust edge crossing helpers ported from S2's `s2edge_crossings`.

use crate::s2::edge_crosser::EdgeCrosser;
use crate::s2::point::{Point, ordered_ccw};

/// Returns a robust cross product direction for unit-length vectors.
pub fn robust_cross_prod(a: &Point, b: &Point) -> Point {
    a.cross(b)
}

#[inline]
fn ref_dir(a: &Point) -> Point {
    a.ortho()
}

/// Returns +1 if AB crosses CD at interior points, 0 for shared vertices
/// between different edges, and -1 otherwise.
pub fn crossing_sign(a: &Point, b: &Point, c: &Point, d: &Point) -> i32 {
    let mut crosser = EdgeCrosser::new_with_c(*a, *b, *c);
    crosser.crossing_sign_with_d(d)
}

/// Returns true if angle ABC contains its vertex B under S2's semi-open rule.
pub fn angle_contains_vertex(a: &Point, b: &Point, c: &Point) -> bool {
    !ordered_ccw(&ref_dir(b), c, a, b)
}

/// Defines whether two edges that share one or more vertices should count as
/// a crossing for point-in-polygon and edge-chain crossing counting.
pub fn vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> bool {
    if a == b || c == d {
        return false;
    }
    if a == c {
        return (b == d) || ordered_ccw(&ref_dir(a), d, b, a);
    }
    if b == d {
        return ordered_ccw(&ref_dir(b), c, a, b);
    }
    if a == d {
        return (b == c) || ordered_ccw(&ref_dir(a), c, b, a);
    }
    if b == c {
        return ordered_ccw(&ref_dir(b), d, a, b);
    }
    false
}

/// Like `vertex_crossing`, but returns crossing orientation:
/// -1 for left-to-right, +1 for right-to-left, 0 otherwise.
pub fn signed_vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> i32 {
    if a == b || c == d {
        return 0;
    }
    if a == c {
        return if (b == d) || ordered_ccw(&ref_dir(a), d, b, a) {
            1
        } else {
            0
        };
    }
    if b == d {
        return if ordered_ccw(&ref_dir(b), c, a, b) { 1 } else { 0 };
    }
    if a == d {
        return if (b == c) || ordered_ccw(&ref_dir(a), c, b, a) {
            -1
        } else {
            0
        };
    }
    if b == c {
        return if ordered_ccw(&ref_dir(b), d, a, b) {
            -1
        } else {
            0
        };
    }
    0
}

/// Convenience wrapper combining `crossing_sign` and `vertex_crossing`.
pub fn edge_or_vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> bool {
    let crossing = crossing_sign(a, b, c, d);
    if crossing < 0 {
        return false;
    }
    if crossing > 0 {
        return true;
    }
    vertex_crossing(a, b, c, d)
}

/// Signed version of `edge_or_vertex_crossing`.
pub fn signed_edge_or_vertex_crossing(a: &Point, b: &Point, c: &Point, d: &Point) -> i32 {
    let crossing = crossing_sign(a, b, c, d);
    if crossing < 0 {
        return 0;
    }
    if crossing > 0 {
        let mut crosser = EdgeCrosser::new_with_c(*a, *b, *c);
        let _ = crosser.crossing_sign_with_d(d);
        return crosser.last_interior_crossing_sign();
    }
    signed_vertex_crossing(a, b, c, d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r3::vector::Vector;
    use crate::s2::point::ORIGIN;

    fn norm(v: Vector) -> Point {
        Point(v).normalize()
    }

    fn test_crossings(a: Point, b: Point, c: Point, d: Point, crossing_sign_want: i32, signed_want: i32) {
        assert_eq!(crossing_sign_want, crossing_sign(&a, &b, &c, &d));
        assert_eq!(signed_want != 0, edge_or_vertex_crossing(&a, &b, &c, &d));
        assert_eq!(signed_want, signed_edge_or_vertex_crossing(&a, &b, &c, &d));
    }

    #[test]
    fn test_crossings_basic() {
        test_crossings(
            norm(Vector::new(1., 2., 1.)),
            norm(Vector::new(1., -3., 0.5)),
            norm(Vector::new(1., -0.5, -3.)),
            norm(Vector::new(0.1, 0.5, 3.)),
            1,
            1,
        );
        test_crossings(
            norm(Vector::new(1., 2., 1.)),
            norm(Vector::new(1., -3., 0.5)),
            norm(Vector::new(-1., 0.5, 3.)),
            norm(Vector::new(-0.1, -0.5, -3.)),
            -1,
            0,
        );
        test_crossings(
            norm(Vector::new(7., -2., 3.)),
            norm(Vector::new(2., 3., 4.)),
            norm(Vector::new(2., 3., 4.)),
            norm(Vector::new(-1., 2., 5.)),
            0,
            -1,
        );
        test_crossings(
            norm(Vector::new(1., 0., 0.)),
            ORIGIN,
            norm(Vector::new(1., -0.1, 1.)),
            norm(Vector::new(1., 1., -0.1)),
            1,
            1,
        );
    }

    #[test]
    fn test_vertex_crossing_symmetry_properties() {
        let a = norm(Vector::new(7., -2., 3.));
        let b = norm(Vector::new(2., 3., 4.));
        let d = norm(Vector::new(-1., 2., 5.));
        assert!(vertex_crossing(&a, &b, &b, &d));
        assert!(vertex_crossing(&b, &a, &b, &d));
        assert!(vertex_crossing(&a, &b, &d, &b));
        assert!(vertex_crossing(&b, &a, &d, &b));
    }

    #[test]
    fn test_angle_contains_vertex() {
        let a = norm(Vector::new(1., 0., 0.));
        let b = norm(Vector::new(0., 1., 0.));
        let c = norm(Vector::new(0., 0., 1.));
        assert!(!angle_contains_vertex(&a, &b, &a));
        assert!(angle_contains_vertex(&ref_dir(&b), &b, &a));
        assert!(!angle_contains_vertex(&a, &b, &ref_dir(&b)));
        assert!(angle_contains_vertex(&a, &b, &c) ^ angle_contains_vertex(&c, &b, &a));
    }
}
