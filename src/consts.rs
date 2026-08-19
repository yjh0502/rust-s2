//! Define the maximum rounding error for arithmetic operations. Depending on the
//! platform the mantissa precision may be different than others, so we choose to
//! use specific values to be consistent across all.
//! The values come from the C++ implementation.

/// EPSILON is a small number that represents a reasonable level of noise between two
/// values that can be considered to be equal.
pub const EPSILON: f64 = 1e-14;

/// DBL_EPSILON is a smaller number for values that require more precision.
pub const DBL_EPSILON: f64 = 2.220446049250313e-16;

#[macro_export]
macro_rules! f64_eq {
    ($x:expr, $y:expr) => {
        ($x - $y).abs() < EPSILON
    };
}

#[macro_export]
macro_rules! assert_f64_eq {
    ($x:expr, $y:expr) => {
        assert!(($x - $y).abs() < EPSILON)
    };
}

#[allow(unused)]
/// f64_eq reports whether the two values are within the default epsilon.
pub fn f64_eq(x: f64, y: f64) -> bool {
    f64_near(x, y, EPSILON)
}

#[allow(unused)]
/// f64_near reports whether the two values are within the specified epsilon.
pub fn f64_near(x: f64, y: f64, eps: f64) -> bool {
    (x - y).abs() <= eps
}

///TODO: to util module?
pub fn remainder(x: f64, y: f64) -> f64 {
    ::libm::remquo(x, y).0
}

/// nextafter returns the next representable f64 after x, moving in the
/// direction of y: y itself if x == y, NaN if either argument is NaN,
/// otherwise the adjacent representable value one step closer to y. This
/// mirrors C's nextafter, built on f64::next_up()/next_down() (stable
/// since Rust 1.86) instead of pulling in a dedicated dependency for it.
pub fn nextafter(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        f64::NAN
    } else if x == y {
        y
    } else if x < y {
        x.next_up()
    } else {
        x.next_down()
    }
}

pub fn clamp<T>(val: T, min: T, max: T) -> T
where
    T: PartialOrd,
{
    if val < min {
        min
    } else if val > max {
        max
    } else {
        val
    }
}

pub fn search_lower_by<F>(len: usize, f: F) -> usize
where
    F: Fn(usize) -> bool,
{
    let mut i = 0;
    let mut j = len;

    while i < j {
        let h = i + (j - i) / 2;
        if !f(h) {
            i = h + 1;
        } else {
            j = h;
        }
    }
    i
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nextafter() {
        // Equal inputs: returns y unchanged, regardless of direction.
        assert_eq!(1.0, nextafter(1.0, 1.0));
        assert_eq!(0.0, nextafter(0.0, 0.0));

        // Moves by exactly one ULP, in the direction of y.
        assert_eq!(1.0f64.next_up(), nextafter(1.0, 2.0));
        assert_eq!(1.0f64.next_down(), nextafter(1.0, 0.0));

        // Crossing zero in either direction steps to the smallest subnormal
        // of the appropriate sign.
        assert_eq!(-5e-324, nextafter(0.0, -1.0));
        assert_eq!(5e-324, nextafter(-0.0, 1.0));

        // NaN in either position propagates to NaN.
        assert!(nextafter(f64::NAN, 1.0).is_nan());
        assert!(nextafter(1.0, f64::NAN).is_nan());
    }
}
