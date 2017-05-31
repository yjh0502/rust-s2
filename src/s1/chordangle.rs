/*
Copyright 2015 Google Inc. All rights reserved.
Copyright 2017 Jihyun Yu. All rights reserved.

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

use std;

use consts::*;
use s1::angle::Angle;

/// ChordAngle represents the angle subtended by a chord (i.e., the straight
/// line segment connecting two points on the sphere). Its representation
/// makes it very efficient for computing and comparing distances, but unlike
/// Angle it is only capable of representing angles between 0 and π radians.
/// Generally, ChordAngle should only be used in loops where many angles need
/// to be calculated and compared. Otherwise it is simpler to use Angle.
///
/// ChordAngle loses some accuracy as the angle approaches π radians.
/// Specifically, the representation of (π - x) radians has an error of about
/// (1e-15 / x), with a maximum error of about 2e-8 radians (about 13cm on the
/// Earth's surface). For comparison, for angles up to π/2 radians (10000km)
/// the worst-case representation error is about 2e-16 radians (1 nanonmeter),
/// which is about the same as Angle.
///
/// ChordAngles are represented by the squared chord length, which can
/// range from 0 to 4. Positive infinity represents an infinite squared length.
#[derive(Debug,PartialEq,PartialOrd,Clone)]
pub struct ChordAngle(pub f64);

/// NEGATIVE represents a chord angle smaller than the zero angle.
/// The only valid operations on a NegativeChordAngle are comparisons and
/// Angle conversions.
pub const NEGATIVE: ChordAngle = ChordAngle(-1f64);

/// RIGHT represents a chord angle of 90 degrees (a "right angle").
pub const RIGHT: ChordAngle = ChordAngle(2f64);

/// STRAIGHT represents a chord angle of 180 degrees (a "straight angle").
/// This is the maximum finite chord angle.
pub const STRAIGHT: ChordAngle = ChordAngle(4f64);

impl<'a> From<&'a Angle> for ChordAngle {
    /// returns a ChordAngle from the given Angle.
    fn from(a: &'a Angle) -> Self {
        if a.0 < 0. {
            NEGATIVE
        } else if a.is_infinite() {
            ChordAngle::inf()
        } else {
            let l = 2. * (0.5 * a.rad().min(std::f64::consts::PI)).sin();
            ChordAngle(l * l)
        }
    }
}
impl From<Angle> for ChordAngle {
    /// returns a ChordAngle from the given Angle.
    fn from(a: Angle) -> Self {
        ChordAngle::from(&a)
    }
}

impl<'a> From<&'a ChordAngle> for Angle {
    /// converts this ChordAngle to an Angle.
    fn from(ca: &'a ChordAngle) -> Self {
        if ca.0 < 0. {
            Angle(-1.)
        } else if ca.is_infinite() {
            Angle::inf()
        } else {
            Angle(2f64 * (0.5 * ca.0.sqrt()).asin())
        }

    }
}
impl From<ChordAngle> for Angle {
    /// converts this ChordAngle to an Angle.
    fn from(ca: ChordAngle) -> Self {
        Angle::from(&ca)
    }
}

impl<'a, 'b> std::ops::Add<&'a ChordAngle> for &'b ChordAngle {
    type Output = ChordAngle;
    /// add adds the other ChordAngle to this one and returns the resulting value.
    /// This method assumes the ChordAngles are not special.
    fn add(self, other: &'a ChordAngle) -> Self::Output {
        // Note that this method (and Sub) is much more efficient than converting
        // the ChordAngle to an Angle and adding those and converting back. It
        // requires only one square root plus a few additions and multiplications.

        if other.0 == 0.0 {
            // Optimization for the common case where b is an error tolerance
            // parameter that happens to be set to zero.
            self.clone()
        } else if self.0 + other.0 >= 4. {
            // Clamp the angle sum to at most 180 degrees.
            STRAIGHT
        } else {
            // Let a and b be the (non-squared) chord lengths, and let c = a+b.
            // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
            // Then the formula below can be derived from c = 2 * sin(A+B) and the
            // relationships   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
            //                 cos(X) = sqrt(1 - sin^2(X))
            let x = self.0 * (1. - 0.25 * other.0);
            let y = other.0 * (1. - 0.25 * self.0);
            ChordAngle(4f64.min(x + y + 2f64 * (x * y).sqrt()))
        }
    }
}

impl std::ops::Add<ChordAngle> for ChordAngle {
    type Output = ChordAngle;
    fn add(self, other: ChordAngle) -> Self::Output {
        &self + &other
    }
}

impl std::ops::Sub<ChordAngle> for ChordAngle {
    type Output = ChordAngle;
    /// sub subtracts the other ChordAngle from this one and returns the resulting
    /// value. This method assumes the ChordAngles are not special.
    fn sub(self, other: ChordAngle) -> Self::Output {
        if other.0 == 0.0 {
            self
        } else if self.0 <= other.0 {
            ChordAngle(0f64)
        } else {
            let x = self.0 * (1. - 0.25 * other.0);
            let y = other.0 * (1. - 0.25 * self.0);
            ChordAngle(0f64.max(x + y - 2. * (x * y).sqrt()))
        }
    }
}

impl ChordAngle {
    /// inf returns a chord angle larger than any finite chord angle.
    /// The only valid operations on an InfChordAngle are comparisons and Angle conversions.
    pub fn inf() -> Self {
        ChordAngle(std::f64::INFINITY)
    }

    /// is_infinite reports whether this ChordAngle is infinite.
    pub fn is_infinite(&self) -> bool {
        self.0.is_infinite()
    }

    /// from_squared_length returns a ChordAngle from the squared chord length.
    /// Note that the argument is automatically clamped to a maximum of 4.0 to
    /// handle possible roundoff errors. The argument must be non-negative.
    pub fn from_squared_length(length2: f64) -> Self {
        if length2 > 4. {
            STRAIGHT
        } else {
            ChordAngle(length2)
        }
    }

    /// expanded returns a new ChordAngle that has been adjusted by the given error
    /// bound (which can be positive or negative). Error should be the value
    /// returned by either MaxPointError or MaxAngleError. For example:
    ///     let a = ChordAngle::from_points(x, y)
    ///     let a1 = a.expanded(a.max_point_error())
    pub fn expanded(&self, e: f64) -> Self {
        // If the angle is special, don't change it. Otherwise clamp it to the valid range.
        if self.is_special() {
            self.clone()
        } else {
            return ChordAngle(0f64.max(4f64.min(self.0 + e)));
        }
    }

    /// is_special reports whether this ChordAngle is one of the special cases.
    pub fn is_special(&self) -> bool {
        self.0 < 0. || self.0.is_infinite()
    }

    /// is_valid reports whether this ChordAngle is valid or not.
    pub fn is_valid(&self) -> bool {
        self.0 >= 0. && self.0 <= 4. || self.is_special()
    }

    /// max_point_error returns the maximum error size for a ChordAngle constructed
    /// from 2 Points x and y, assuming that x and y are normalized to within the
    /// bounds guaranteed by s2.Point.Normalize. The error is defined with respect to
    /// the true distance after the points are projected to lie exactly on the sphere.
    pub fn max_point_error(&self) -> f64 {
        // There is a relative error of (2.5*DBL_EPSILON) when computing the squared
        // distance, plus an absolute error of (16 * DBL_EPSILON**2) because the
        // lengths of the input points may differ from 1 by up to (2*DBL_EPSILON) each.
        2.5 * DBL_EPSILON * self.0 + 16. * DBL_EPSILON * DBL_EPSILON
    }

    /// max_angle_error returns the maximum error for a ChordAngle constructed
    /// as an Angle distance.
    pub fn max_angle_error(&self) -> f64 {
        DBL_EPSILON * self.0
    }

    /// sin returns the sine of this chord angle. This method is more efficient
    /// than converting to Angle and performing the computation.
    pub fn sin(&self) -> f64 {
        self.sin2().sqrt()
    }

    /// sin2 returns the square of the sine of this chord angle.
    /// It is more efficient than Sin.
    pub fn sin2(&self) -> f64 {
        // Let a be the (non-squared) chord length, and let A be the corresponding
        // half-angle (a = 2*sin(A)).  The formula below can be derived from:
        //   sin(2*A) = 2 * sin(A) * cos(A)
        //   cos^2(A) = 1 - sin^2(A)
        // This is much faster than converting to an angle and computing its sine.
        self.0 * (1. - 0.25 * self.0)
    }

    /// cos returns the cosine of this chord angle. This method is more efficient
    /// than converting to Angle and performing the computation.
    pub fn cos(&self) -> f64 {
        // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
        1.0 - 0.5 * self.0
    }

    /// tan returns the tangent of this chord angle.
    pub fn tan(&self) -> f64 {
        self.sin() / self.cos()
    }
}

/*
package s1

import (
	"math"
	"testing"
)

func TestChordAngleBasics(t *testing.T) {
	var zeroChord ChordAngle
	tests := []struct {
		a, b     ChordAngle
		lessThan bool
		equal    bool
	}{
		{NegativeChordAngle, NegativeChordAngle, false, true},
		{NegativeChordAngle, zeroChord, true, false},
		{NegativeChordAngle, StraightChordAngle, true, false},
		{NegativeChordAngle, InfChordAngle(), true, false},

		{zeroChord, zeroChord, false, true},
		{zeroChord, StraightChordAngle, true, false},
		{zeroChord, InfChordAngle(), true, false},

		{StraightChordAngle, StraightChordAngle, false, true},
		{StraightChordAngle, InfChordAngle(), true, false},

		{InfChordAngle(), InfChordAngle(), false, true},
		{InfChordAngle(), InfChordAngle(), false, true},
	}

	for _, test := range tests {
		if got := test.a < test.b; got != test.lessThan {
			t.Errorf("%v should be less than %v", test.a, test.b)
		}
		if got := test.a == test.b; got != test.equal {
			t.Errorf("%v should be equal to %v", test.a, test.b)
		}
	}
}

func TestChordAngleIsFunctions(t *testing.T) {
	var zeroChord ChordAngle
	tests := []struct {
		have       ChordAngle
		isNegative bool
		isZero     bool
		isInf      bool
		isSpecial  bool
	}{
		{zeroChord, false, true, false, false},
		{NegativeChordAngle, true, false, false, true},
		{zeroChord, false, true, false, false},
		{StraightChordAngle, false, false, false, false},
		{InfChordAngle(), false, false, true, true},
	}

	for _, test := range tests {
		if got := test.have < 0; got != test.isNegative {
			t.Errorf("%v.isNegative() = %t, want %t", test.have, got, test.isNegative)
		}
		if got := test.have == 0; got != test.isZero {
			t.Errorf("%v.isZero() = %t, want %t", test.have, got, test.isZero)
		}
		if got := test.have.isInf(); got != test.isInf {
			t.Errorf("%v.isInf() = %t, want %t", test.have, got, test.isInf)
		}
		if got := test.have.isSpecial(); got != test.isSpecial {
			t.Errorf("%v.isSpecial() = %t, want %t", test.have, got, test.isSpecial)
		}
	}
}

func TestChordAngleFromAngle(t *testing.T) {
	for _, angle := range []float64{0, 1, -1, math.Pi} {
		if got := ChordAngleFromAngle(Angle(angle)).Angle().Radians(); got != angle {
			t.Errorf("ChordAngleFromAngle(Angle(%v)) = %v, want %v", angle, got, angle)
		}
	}

	if got := ChordAngleFromAngle(Angle(math.Pi)); got != StraightChordAngle {
		t.Errorf("a ChordAngle from an Angle of π = %v, want %v", got, StraightChordAngle)
	}

	if InfAngle() != ChordAngleFromAngle(InfAngle()).Angle() {
		t.Errorf("converting infinite Angle to ChordAngle should yield infinite Angle")
	}
}

func TestChordAngleArithmetic(t *testing.T) {
	var (
		zero      ChordAngle
		degree30  = ChordAngleFromAngle(30 * Degree)
		degree60  = ChordAngleFromAngle(60 * Degree)
		degree90  = ChordAngleFromAngle(90 * Degree)
		degree120 = ChordAngleFromAngle(120 * Degree)
		degree180 = StraightChordAngle
	)

	addTests := []struct {
		a, b ChordAngle
		want ChordAngle
	}{
		{zero, zero, zero},
		{degree60, zero, degree60},
		{zero, degree60, degree60},
		{degree30, degree60, degree90},
		{degree60, degree30, degree90},
		{degree180, zero, degree180},
		{degree60, degree30, degree90},
		{degree90, degree90, degree180},
		{degree120, degree90, degree180},
		{degree120, degree120, degree180},
		{degree30, degree180, degree180},
		{degree180, degree180, degree180},
	}

	subTests := []struct {
		a, b ChordAngle
		want ChordAngle
	}{
		{zero, zero, zero},
		{degree60, degree60, zero},
		{degree180, degree180, zero},
		{zero, degree60, zero},
		{degree30, degree90, zero},
		{degree90, degree30, degree60},
		{degree90, degree60, degree30},
		{degree180, zero, degree180},
	}

	for _, test := range addTests {
		if got := float64(test.a.Add(test.b)); !float64Eq(got, float64(test.want)) {
			t.Errorf("%v.Add(%v) = %0.24f, want %0.24f", test.a.Angle().Degrees(), test.b.Angle().Degrees(), got, test.want)
		}
	}
	for _, test := range subTests {
		if got := float64(test.a.Sub(test.b)); !float64Eq(got, float64(test.want)) {
			t.Errorf("%v.Sub(%v) = %0.24f, want %0.24f", test.a.Angle().Degrees(), test.b.Angle().Degrees(), got, test.want)
		}
	}
}

func TestChordAngleTrigonometry(t *testing.T) {
	// Because of the way the math works out, the 9/10th's case has slightly more
	// difference than all the other computations, so this gets a more generous
	// epsilon to deal with that.
	const epsilon = 1e-14
	const iters = 40
	for iter := 0; iter <= iters; iter++ {
		radians := math.Pi * float64(iter) / float64(iters)
		angle := ChordAngleFromAngle(Angle(radians))
		if !float64Near(math.Sin(radians), angle.Sin(), epsilon) {
			t.Errorf("(%d/%d)*π. %v.Sin() = %v, want %v", iter, iters, angle, angle.Sin(), math.Sin(radians))
		}
		if !float64Near(math.Cos(radians), angle.Cos(), epsilon) {
			t.Errorf("(%d/%d)*π. %v.Cos() = %v, want %v", iter, iters, angle, angle.Cos(), math.Cos(radians))
		}
		// Since tan(x) is unbounded near pi/4, we map the result back to an
		// angle before comparing. The assertion is that the result is equal to
		// the tangent of a nearby angle.
		if !float64Near(math.Atan(math.Tan(radians)), math.Atan(angle.Tan()), 1e-14) {
			t.Errorf("(%d/%d)*π. %v.Tan() = %v, want %v", iter, iters, angle, angle.Tan(), math.Tan(radians))
		}
	}

	// Unlike Angle, ChordAngle can represent 90 and 180 degrees exactly.
	angle90 := ChordAngleFromSquaredLength(2)
	angle180 := ChordAngleFromSquaredLength(4)
	if !float64Eq(1, angle90.Sin()) {
		t.Errorf("%v.Sin() = %v, want 1", angle90, angle90.Sin())
	}
	if !float64Eq(0, angle90.Cos()) {
		t.Errorf("%v.Cos() = %v, want 0", angle90, angle90.Cos())
	}
	if !math.IsInf(angle90.Tan(), 0) {
		t.Errorf("%v.Tan() should be infinite, but was not.", angle90)
	}
	if !float64Eq(0, angle180.Sin()) {
		t.Errorf("%v.Sin() = %v, want 0", angle180, angle180.Sin())
	}
	if !float64Eq(-1, angle180.Cos()) {
		t.Errorf("%v.Cos() = %v, want -1", angle180, angle180.Cos())
	}
	if !float64Eq(0, angle180.Tan()) {
		t.Errorf("%v.Tan() = %v, want 0", angle180, angle180.Tan())
	}
}

func TestChordAngleExpanded(t *testing.T) {
	var zero ChordAngle

	tests := []struct {
		have ChordAngle
		add  float64
		want ChordAngle
	}{
		{NegativeChordAngle, 5, NegativeChordAngle.Expanded(5)},
		{InfChordAngle(), -5, InfChordAngle()},
		{StraightChordAngle, 5, ChordAngleFromSquaredLength(5)},
		{zero, -5, zero},
		{ChordAngleFromSquaredLength(1.25), 0.25, ChordAngleFromSquaredLength(1.5)},
		{ChordAngleFromSquaredLength(0.75), 0.25, ChordAngleFromSquaredLength(1)},
	}

	for _, test := range tests {
		if got := test.have.Expanded(test.add); got != test.want {
			t.Errorf("%v.Expanded(%v) = %v, want %v", test.have, test.add, got, test.want)
		}
	}
}
*/
