/*
Copyright 2026 The rust-s2 Authors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
*/

//! Stateful edge crosser ported from S2's `s2edge_crosser`.

use crate::consts::DBL_EPSILON;
use crate::r3::vector::Vector;
use crate::s2::edge_crossings::{signed_vertex_crossing, vertex_crossing};
use crate::s2::point::Point;
use crate::s2::predicates::{Direction, robust_sign, sign, triage_sign};

#[inline]
fn direction_to_int(d: Direction) -> i32 {
    match d {
        Direction::Clockwise => -1,
        Direction::Indeterminate => 0,
        Direction::CounterClockwise => 1,
    }
}

#[inline]
fn triage_sign_with_cross(cross: &Vector, c: &Point) -> i32 {
    const MAX_DETERMINANT_ERROR: f64 = 1.8274 * DBL_EPSILON;
    let det = cross.dot(&c.0);
    if det > MAX_DETERMINANT_ERROR {
        1
    } else if det < -MAX_DETERMINANT_ERROR {
        -1
    } else {
        0
    }
}

#[inline]
fn robust_sign_int(a: &Point, b: &Point, c: &Point) -> i32 {
    if a == b || b == c || c == a {
        return 0;
    }
    match robust_sign(a, b, c) {
        Direction::CounterClockwise => 1,
        Direction::Clockwise => -1,
        // exact_sign isn't implemented in this crate yet; use a deterministic
        // fallback that preserves non-zero orientation for distinct points.
        Direction::Indeterminate => {
            if sign(a, b, c) {
                1
            } else {
                -1
            }
        }
    }
}

#[inline]
fn sign_with_cross(a: &Point, b: &Point, c: &Point, cross: &Vector) -> i32 {
    if a == b || b == c || c == a {
        return 0;
    }
    let triage = triage_sign_with_cross(cross, c);
    if triage != 0 {
        triage
    } else {
        robust_sign_int(a, b, c)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct EdgeCrosser {
    a: Point,
    b: Point,
    a_cross_b: Vector,
    have_tangents: bool,
    a_tangent: Point,
    b_tangent: Point,
    c: Option<Point>,
    acb: i32,
    bda: i32,
}

impl EdgeCrosser {
    pub fn new(a: Point, b: Point) -> Self {
        Self {
            a,
            b,
            a_cross_b: a.0.cross(&b.0),
            have_tangents: false,
            a_tangent: Point::default(),
            b_tangent: Point::default(),
            c: None,
            acb: 0,
            bda: 0,
        }
    }

    pub fn new_with_c(a: Point, b: Point, c: Point) -> Self {
        let mut out = Self::new(a, b);
        out.restart_at(c);
        out
    }

    pub fn a(&self) -> Point {
        self.a
    }

    pub fn b(&self) -> Point {
        self.b
    }

    pub fn c(&self) -> Option<Point> {
        self.c
    }

    pub fn init(&mut self, a: Point, b: Point) {
        self.a = a;
        self.b = b;
        self.a_cross_b = a.0.cross(&b.0);
        self.have_tangents = false;
        self.c = None;
        self.acb = 0;
        self.bda = 0;
    }

    pub fn restart_at(&mut self, c: Point) {
        self.c = Some(c);
        self.acb = -direction_to_int(triage_sign(&self.a, &self.b, &c));
    }

    pub fn crossing_sign(&mut self, c: &Point, d: &Point) -> i32 {
        if self.c != Some(*c) {
            self.restart_at(*c);
        }
        self.crossing_sign_with_d(d)
    }

    pub fn crossing_sign_with_d(&mut self, d: &Point) -> i32 {
        let bda = direction_to_int(triage_sign(&self.a, &self.b, d));
        if self.acb == -bda && bda != 0 {
            self.c = Some(*d);
            self.acb = -bda;
            return -1;
        }
        self.bda = bda;
        self.crossing_sign_internal(*d)
    }

    pub fn edge_or_vertex_crossing(&mut self, c: &Point, d: &Point) -> bool {
        if self.c != Some(*c) {
            self.restart_at(*c);
        }
        self.edge_or_vertex_crossing_with_d(d)
    }

    pub fn edge_or_vertex_crossing_with_d(&mut self, d: &Point) -> bool {
        let c = self.c.expect("restart_at must be called before chain crossing");
        let crossing = self.crossing_sign_with_d(d);
        if crossing < 0 {
            return false;
        }
        if crossing > 0 {
            return true;
        }
        vertex_crossing(&self.a, &self.b, &c, d)
    }

    pub fn signed_edge_or_vertex_crossing(&mut self, c: &Point, d: &Point) -> i32 {
        if self.c != Some(*c) {
            self.restart_at(*c);
        }
        self.signed_edge_or_vertex_crossing_with_d(d)
    }

    pub fn signed_edge_or_vertex_crossing_with_d(&mut self, d: &Point) -> i32 {
        let c = self.c.expect("restart_at must be called before chain crossing");
        let crossing = self.crossing_sign_with_d(d);
        if crossing < 0 {
            return 0;
        }
        if crossing > 0 {
            return self.last_interior_crossing_sign();
        }
        signed_vertex_crossing(&self.a, &self.b, &c, d)
    }

    pub fn last_interior_crossing_sign(&self) -> i32 {
        self.acb
    }

    fn crossing_sign_internal(&mut self, d: Point) -> i32 {
        let result = self.crossing_sign_internal2(&d);
        self.c = Some(d);
        self.acb = -self.bda;
        result
    }

    fn crossing_sign_internal2(&mut self, d: &Point) -> i32 {
        let c = self.c.expect("restart_at must be called before chain crossing");

        if !self.have_tangents {
            let norm = self.a.cross(&self.b);
            self.a_tangent = Point(self.a.0.cross(&norm.0));
            self.b_tangent = Point(norm.0.cross(&self.b.0));
            self.have_tangents = true;
        }

        let k_error = (1.5 + 1.0 / (3.0f64).sqrt()) * DBL_EPSILON;
        if (c.0.dot(&self.a_tangent.0) > k_error && d.0.dot(&self.a_tangent.0) > k_error)
            || (c.0.dot(&self.b_tangent.0) > k_error && d.0.dot(&self.b_tangent.0) > k_error)
        {
            return -1;
        }

        if self.a == c || self.a == *d || self.b == c || self.b == *d {
            return 0;
        }
        if self.a == self.b || c == *d {
            return -1;
        }

        if self.acb == 0 {
            self.acb = -robust_sign_int(&self.a, &self.b, &c);
        }
        if self.bda == 0 {
            self.bda = robust_sign_int(&self.a, &self.b, d);
        }
        if self.bda != self.acb {
            return -1;
        }

        let c_cross_d = c.0.cross(&d.0);
        let cbd = -sign_with_cross(&c, d, &self.b, &c_cross_d);
        if cbd != self.acb {
            return -1;
        }
        let dac = sign_with_cross(&c, d, &self.a, &c_cross_d);
        if dac != self.acb { -1 } else { 1 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r3::vector::Vector;
    use crate::s2::edge_crossings::{crossing_sign, edge_or_vertex_crossing};

    fn norm(v: Vector) -> Point {
        Point(v).normalize()
    }

    #[test]
    fn test_crossing_sign_matches_stateless() {
        let a = norm(Vector::new(1., 2., 1.));
        let b = norm(Vector::new(1., -3., 0.5));
        let c = norm(Vector::new(1., -0.5, -3.));
        let d = norm(Vector::new(0.1, 0.5, 3.));

        let mut crosser = EdgeCrosser::new_with_c(a, b, c);
        let got = crosser.crossing_sign_with_d(&d);
        assert_eq!(crossing_sign(&a, &b, &c, &d), got);
        assert_eq!(1, got);
    }

    #[test]
    fn test_edge_chain_methods() {
        let a = norm(Vector::new(1., 2., 1.));
        let b = norm(Vector::new(1., -3., 0.5));
        let chain = [
            norm(Vector::new(1., -0.5, -3.)),
            norm(Vector::new(0.1, 0.5, 3.)),
            norm(Vector::new(2., 0.5, -3.)),
        ];

        let mut crosser = EdgeCrosser::new_with_c(a, b, chain[0]);
        let mut c = chain[0];
        for d in chain.iter().skip(1) {
            let got = crosser.edge_or_vertex_crossing_with_d(d);
            let want = edge_or_vertex_crossing(&a, &b, &c, d);
            assert_eq!(want, got);
            c = *d;
        }
    }

    #[test]
    fn test_reuse_via_init() {
        let a = norm(Vector::new(1., 0., 0.));
        let b = norm(Vector::new(0., 1., 0.));
        let c = norm(Vector::new(0., 0., 1.));
        let d = norm(Vector::new(1., 1., 1.));

        let mut crosser = EdgeCrosser::new(a, b);
        crosser.restart_at(c);
        let before = crosser.crossing_sign_with_d(&d);

        crosser.init(c, d);
        crosser.restart_at(a);
        let after = crosser.crossing_sign_with_d(&b);
        assert_eq!(before, after);
    }
}
