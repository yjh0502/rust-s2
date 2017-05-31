use std;

use consts::EPSILON;
use s1::angle::Angle;

#[derive(Clone,PartialEq,PartialOrd,Debug)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl std::ops::Add<Vector> for Vector {
    type Output = Vector;
    fn add(self, other: Vector) -> Self::Output {
        &self + &other
    }
}
impl<'a, 'b> std::ops::Add<&'b Vector> for &'a Vector {
    type Output = Vector;
    fn add(self, other: &'b Vector) -> Self::Output {
        Vector {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl std::ops::Sub<Vector> for Vector {
    type Output = Vector;
    fn sub(self, other: Vector) -> Self::Output {
        &self - &other
    }
}
impl<'a, 'b> std::ops::Sub<&'b Vector> for &'a Vector {
    type Output = Vector;
    fn sub(self, other: &'b Vector) -> Self::Output {
        Vector {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl std::ops::Mul<Vector> for Vector {
    type Output = Vector;
    fn mul(self, other: Vector) -> Self::Output {
        Vector {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z,
        }
    }
}

impl<'a> std::ops::Mul<f64> for &'a Vector {
    type Output = Vector;
    fn mul(self, m: f64) -> Self::Output {
        Vector {
            x: self.x * m,
            y: self.y * m,
            z: self.z * m,
        }
    }
}
impl std::ops::Mul<f64> for Vector {
    type Output = Vector;
    fn mul(self, m: f64) -> Self::Output {
        &self * m
    }
}

impl Vector {
    pub fn xyz(x: f64, y: f64, z: f64) -> Self {
        Vector { x: x, y: y, z: z }
    }

    pub fn approx_eq(&self, other: &Vector) -> bool {
        (self.x - other.x).abs() < EPSILON && (self.y - other.y).abs() < EPSILON && (self.z - other.z).abs() < EPSILON
    }

    pub fn norm(&self) -> f64 {
        self.norm2().sqrt()
    }

    pub fn norm2(&self) -> f64 {
        self.dot(self)
    }

    pub fn normalize(&self) -> Self {
        if self.x == 0. && self.y == 0. && self.z == 0. {
            self.clone()
        } else {
            self.clone() * (1.0 / self.norm())
        }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn is_unit(&self) -> bool {
        const EPSILON2: f64 = 5e-14;
        (self.norm2() - 1.).abs() < EPSILON2
    }

    pub fn abs(&self) -> Self {
        Vector {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }

    pub fn cross(&self, other: &Self) -> Self {
        Vector {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn distance(&self, other: &Self) -> f64 {
        (self - other).norm()
    }

    pub fn angle(&self, other: &Self) -> Angle {
        Angle(self.cross(other).norm().atan2(self.dot(other)))
    }

    pub fn ortho(&self) -> Self {
        let mut ov = Self {
            x: 0.012,
            y: 0.0053,
            z: 0.00457,
        };
        match self.largest_component() {
            Axis::X => ov.z = 1.0,
            Axis::Y => ov.x = 1.0,
            Axis::Z => ov.y = 1.0,
        };
        self.cross(&ov).normalize()
    }

    pub fn largest_component(&self) -> Axis {
        let a = self.abs();
        if a.x > a.y {
            if a.x > a.z { Axis::X } else { Axis::Z }
        } else {
            if a.y > a.z { Axis::Y } else { Axis::Z }
        }
    }

    // smallest_component returns the axis that represents the smallest component in this vector.
    pub fn smallest_component(&self) -> Axis {
        let t = self.abs();
        if t.x < t.y {
            if t.x < t.z { Axis::X } else { Axis::Z }
        } else {
            if t.y < t.z { Axis::Y } else { Axis::Z }
        }
    }
}

#[derive(PartialEq,Eq,Debug)]
pub enum Axis {
    X,
    Y,
    Z,
}


#[cfg(test)]
mod tests {
    use super::*;
    use consts::*;
    use std::cmp::Ordering;
    use std::f64::consts::PI;

    macro_rules! V {
        ($x: expr, $y: expr, $z: expr) => {
            Vector{x: $x, y: $y , z: $z }
        }
    }

    #[test]
    fn test_vector_norm() {
        assert_eq!(V!(0., 0., 0.).norm(), 0.);
        assert_eq!(V!(0., 1., 0.).norm(), 1.);
        assert_eq!(V!(3., -4., 12.).norm(), 13.);
        assert_eq!(V!(1., 1e-16, 1e-32).norm(), 1.);
    }

    #[test]
    fn test_vector_norm2() {
        assert_eq!(V!(0., 0., 0.).norm2(), 0.);
        assert_eq!(V!(0., 1., 0.).norm2(), 1.);
        assert_eq!(V!(1., 1., 1.).norm2(), 3.);
        assert_eq!(V!(1., 2., 3.).norm2(), 14.);
        assert_eq!(V!(3., -4., 12.).norm2(), 169.);
        assert_eq!(V!(1., 1e-16, 1e-32).norm2(), 1.);
    }

    fn test_vec_norm(v: Vector) {
        let n = v.normalize();
        assert!(f64_eq(v.x * n.y, v.y * n.x));
        assert!(f64_eq(v.x * n.z, v.z * n.x));

        assert!(f64_eq(n.norm(), 1.));
    }

    #[test]
    fn test_vector_normalize() {
        test_vec_norm(V!(1., 0., 0.));
        test_vec_norm(V!(0., 1., 0.));
        test_vec_norm(V!(0., 0., 1.));
        test_vec_norm(V!(1., 1., 1.));
        test_vec_norm(V!(1., 1e-16, 1e-32));
        test_vec_norm(V!(12.34, 56.78, 91.01));
    }

    #[test]
    fn test_vector_is_unit() {
        assert_eq!(false, V!(0., 0., 0.).is_unit());
        assert_eq!(true, V!(0., 1., 0.).is_unit());
        assert_eq!(true, V!(1. + 2. * EPSILON, 0., 0.).is_unit());
        assert_eq!(true, V!(1. * (1. + EPSILON), 0., 0.).is_unit());
        assert_eq!(false, V!(1., 1., 1.).is_unit());
        assert_eq!(true, V!(1., 1e-16, 1e-32).is_unit());
    }

    fn test_vector_dot_case(expected: f64, v1: Vector, v2: Vector) {
        assert!(f64_eq(expected, v1.dot(&v2)));
        assert!(f64_eq(expected, v2.dot(&v1)));
    }

    #[test]
    fn test_vector_dot() {
        test_vector_dot_case(1., V!(1., 0., 0.), V!(1., 0., 0.));
        test_vector_dot_case(0., V!(1., 0., 0.), V!(0., 1., 0.));
        test_vector_dot_case(0., V!(1., 0., 0.), V!(0., 1., 1.));
        test_vector_dot_case(-3., V!(1., 1., 1.), V!(-1., -1., -1.));
        test_vector_dot_case(-1.9, V!(1., 2., 2.), V!(-0.3, 0.4, -1.2));
    }

    #[test]
    fn test_vector_cross() {
        assert!(V!(1., 0., 0.)
                    .cross(&V!(1., 0., 0.))
                    .approx_eq(&V!(0., 0., 0.)));
        assert!(V!(1., 0., 0.)
                    .cross(&V!(0., 1., 0.))
                    .approx_eq(&V!(0., 0., 1.)));
        assert!(V!(0., 1., 0.)
                    .cross(&V!(1., 0., 0.))
                    .approx_eq(&V!(0., 0., -1.)));
        assert!(V!(1., 2., 3.)
                    .cross(&V!(-4., 5., -6.))
                    .approx_eq(&V!(-27., -6., 13.)));
    }

    #[test]
    fn test_vector_add() {
        assert!((V!(0., 0., 0.) + V!(0., 0., 0.)).approx_eq(&V!(0., 0., 0.)));
        assert!((V!(1., 0., 0.) + V!(0., 0., 0.)).approx_eq(&V!(1., 0., 0.)));
        assert!((V!(1., 2., 3.) + V!(4., 5., 7.)).approx_eq(&V!(5., 7., 10.)));
        assert!((V!(1., -3., 5.) + V!(1., -6., -6.)).approx_eq(&V!(2., -9., -1.)));
    }

    #[test]
    fn test_vector_sub() {
        assert!((V!(0., 0., 0.) - V!(0., 0., 0.)).approx_eq(&V!(0., 0., 0.)));
        assert!((V!(1., 0., 0.) - V!(0., 0., 0.)).approx_eq(&V!(1., 0., 0.)));
        assert!((V!(1., 2., 3.) - V!(4., 5., 7.)).approx_eq(&V!(-3., -3., -4.)));
        assert!((V!(1., -3., 5.) - V!(1., -6., -6.)).approx_eq(&V!(0., 3., 11.)));
    }

    #[test]
    fn test_vector_distance() {
        assert!(f64_eq(V!(1., 0., 0.).distance(&V!(1., 0., 0.)), 0.));
        assert!(f64_eq(V!(1., 0., 0.).distance(&V!(0., 1., 0.)), 1.41421356237310));
        assert!(f64_eq(V!(1., 0., 0.).distance(&V!(0., 1., 1.)), 1.73205080756888));
        assert!(f64_eq(V!(1., 1., 1.).distance(&V!(-1., -1., -1.)),
                       3.46410161513775));
        assert!(f64_eq(V!(1., 2., 2.).distance(&V!(-0.3, 0.4, -1.2)),
                       3.80657326213486));
    }

    #[test]
    fn test_vector_mul() {
        assert!((V!(0., 0., 0.) * 3.).approx_eq(&V!(0., 0., 0.)));
        assert!((V!(1., 0., 0.) * 1.).approx_eq(&V!(1., 0., 0.)));
        assert!((V!(1., 0., 0.) * 0.).approx_eq(&V!(0., 0., 0.)));
        assert!((V!(1., 0., 0.) * 3.).approx_eq(&V!(3., 0., 0.)));
        assert!((V!(1., -3., 5.) * -1.).approx_eq(&V!(-1., 3., -5.)));
        assert!((V!(1., -3., 5.) * 2.).approx_eq(&V!(2., -6., 10.)));
    }

    #[test]
    fn test_vector_angle() {
        assert!(f64_eq(V!(1., 0., 0.).angle(&V!(1., 0., 0.)).0, 0.));
        assert!(f64_eq(V!(1., 0., 0.).angle(&V!(0., 1., 0.)).0, PI / 2.));
        assert!(f64_eq(V!(1., 0., 0.).angle(&V!(0., 1., 1.)).0, PI / 2.));
        assert!(f64_eq(V!(1., 0., 0.).angle(&V!(-1., 0., 0.)).0, PI));
        assert!(f64_eq(V!(1., 2., 3.).angle(&V!(2., 3., -1.)).0, 1.2055891055045298));
    }

    fn test_vector_ortho_case(v: Vector) {
        assert!(f64_eq(v.dot(&v.ortho()), 0.));
        assert!(f64_eq(v.ortho().norm(), 1.));
    }

    #[test]
    fn test_vector_ortho() {
        test_vector_ortho_case(V!(1., 0., 0.));
        test_vector_ortho_case(V!(1., 1., 0.));
        test_vector_ortho_case(V!(1., 2., 3.));
        test_vector_ortho_case(V!(1., -2., -5.));
        test_vector_ortho_case(V!(0.012, 0.0053, 0.00457));
        test_vector_ortho_case(V!(-0.012, -1., -0.00457));
    }

    fn test_vector_identities_case(v1: Vector, v2: Vector) {
        let a1 = v1.angle(&v2).0;
        let a2 = v2.angle(&v1).0;
        let c1 = v1.cross(&v2);
        let c2 = v2.cross(&v1);
        let d1 = v1.dot(&v2);
        let d2 = v2.dot(&v1);

        // angle commuts
        assert!(f64_eq(a1, a2));
        // dot commutes
        assert!(f64_eq(d1, d2));
        // cross anti-commuts
        assert!(c1.approx_eq(&(c2.clone() * -1.)));

        // cross is orthogonal to original vectors
        assert!(f64_eq(v1.dot(&c1), 0.));
        assert!(f64_eq(v2.dot(&c1), 0.));
        assert!(f64_eq(v1.dot(&c2), 0.));
        assert!(f64_eq(v2.dot(&c2), 0.));
    }

    #[test]
    fn test_vector_identities() {
        test_vector_identities_case(V!(0., 0., 0.), V!(0., 0., 0.));
        test_vector_identities_case(V!(0., 0., 0.), V!(0., 1., 2.));
        test_vector_identities_case(V!(1., 0., 0.), V!(0., 1., 0.));
        test_vector_identities_case(V!(1., 0., 0.), V!(0., 1., 1.));
        test_vector_identities_case(V!(1., 1., 1.), V!(-1., -1., -1.));
        test_vector_identities_case(V!(1., 2., 2.), V!(-0.3, 0.4, -1.2));
    }

    fn test_ls(v: Vector, largest: Axis, smallest: Axis) {
        assert_eq!(v.largest_component(), largest);
        assert_eq!(v.smallest_component(), smallest);
    }

    #[test]
    fn test_vector_largest_smallest_components() {
        test_ls(V!(0., 0., 0.), Axis::Z, Axis::Z);
        test_ls(V!(1., 0., 0.), Axis::X, Axis::Z);
        test_ls(V!(1., -1., 0.), Axis::Y, Axis::Z);
        test_ls(V!(-1., -1.1, -1.1), Axis::Z, Axis::X);
        test_ls(V!(0.5, -0.4, -0.5), Axis::Z, Axis::Y);
        test_ls(V!(1e-15, 1e-14, 1e-13), Axis::Z, Axis::X);
    }

    fn test_cmp(v1: Vector, v2: Vector, expected: Ordering) {
        assert_eq!(v1.partial_cmp(&v2), Some(expected));
    }

    #[test]
    fn test_vector_cmp() {
        // let's hope derived PartialCmp compares element in order
        test_cmp(V!(0., 0., 0.), V!(0., 0., 0.), Ordering::Equal);
        test_cmp(V!(0., 0., 0.), V!(1., 0., 0.), Ordering::Less);
        test_cmp(V!(0., 1., 0.), V!(0., 0., 0.), Ordering::Greater);

        test_cmp(V!(1., 2., 3.), V!(3., 2., 1.), Ordering::Less);
        test_cmp(V!(-1., 0., 0.), V!(0., 0., -1.), Ordering::Less);
        test_cmp(V!(8., 6., 4.), V!(7., 5., 3.), Ordering::Greater);
        test_cmp(V!(-1., -0.5, 0.), V!(0., 0., 0.1), Ordering::Less);
        test_cmp(V!(1., 2., 3.), V!(2., 3., 4.), Ordering::Less);
        test_cmp(V!(1.23, 4.56, 7.89), V!(1.23, 4.56, 7.89), Ordering::Equal);
    }
}
