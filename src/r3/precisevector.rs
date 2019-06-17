/*
Copyright 2019 Alexander Haynes. All rights reserved.

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

extern crate bigdecimal;

//use std;
use r3;
use std::str::FromStr;
use r3::precisevector::bigdecimal::ToPrimitive;

pub fn prec_str(s: String) -> bigdecimal::BigDecimal {
	let f = bigdecimal::BigDecimal::from_str(&s).unwrap();
	return f
}

pub fn prec_int(i: i64) -> bigdecimal::BigDecimal {
	return bigdecimal::BigDecimal::from(i);
}

pub fn prec_float(f: f64) -> bigdecimal::BigDecimal {
	return bigdecimal::BigDecimal::from(f);
}


// PreciseVector represents a point in ℝ³ using high-precision values.
// Note that this is NOT a complete implementation because there are some
// operations that Vector supports that are not feasible with arbitrary precision
// math. (e.g., methods that need divison like Normalize, or methods needing a
// square root operation such as Norm)
#[derive(Clone)]
pub struct PreciseVector {
	x: bigdecimal::BigDecimal, 
	y: bigdecimal::BigDecimal, 
	z: bigdecimal::BigDecimal
}

// PreciseVectorFromVector creates a high precision vector from the given Vector.
pub fn precise_vector_from_vector(v: r3::vector::Vector) -> PreciseVector {
	return new_precise_vector(v.x, v.y, v.z)
}

// NewPreciseVector creates a high precision vector from the given floating point values.
pub fn new_precise_vector(x: f64, y: f64, z: f64) -> PreciseVector {
	return PreciseVector{
		x: prec_float(x),
		y: prec_float(y),
		z: prec_float(z),
	}
}

impl PreciseVector {

	// Vector returns this precise vector converted to a Vector.
	pub fn vector(&self) -> r3::vector::Vector {
		// The accuracy flag is ignored on these conversions back to float64.
		let x = self.x.to_f64().unwrap();
		let y = self.y.to_f64().unwrap();
		let z = self.z.to_f64().unwrap();
		return r3::vector::Vector{x, y, z}.normalize()
	}

	// Equal reports whether v and ov are equal.
	pub fn equal(&self, ov: PreciseVector) -> bool {
		return self.x == ov.x && self.y == ov.y && self.z == ov.z
	}

	// Norm2 returns the square of the norm.
	pub fn norm2(&self) -> bigdecimal::BigDecimal {
		return self.dot(&self.clone());
	}

	// IsUnit reports whether this vector is of unit length.
	pub fn is_unit(&self) -> bool {
		return self.norm2() == prec_int(1_i64);
	}

	// Abs returns the vector with nonnegative components.
	pub fn abs(&self) -> PreciseVector {
		return PreciseVector{
			x: bigdecimal::BigDecimal::abs(&self.x),
			y: bigdecimal::BigDecimal::abs(&self.y),
			z: bigdecimal::BigDecimal::abs(&self.z)
		}
	}

	// Add returns the standard vector sum of v and ov.
	pub fn add(&self, ov: PreciseVector) -> PreciseVector {
		let v = self.clone();
		return PreciseVector{
			x: v.x +ov.x,
			y: v.y +ov.y,
			z: v.z +ov.z
		}
	}

	// Sub returns the standard vector difference of v and ov.
	pub fn sub(&self, ov: PreciseVector) -> PreciseVector {
		let v = self.clone();
		return PreciseVector{
			x: v.x - ov.x,
			y: v.y - ov.y,
			z: v.z - ov.z,
		}
	}

	// Mul returns the standard scalar product of v and f.
	pub fn mul(&self, f: bigdecimal::BigDecimal) -> PreciseVector {
		let v = self.clone();
		return PreciseVector{
			x: v.x * f.clone(),
			y: v.y * f.clone(),
			z: v.z * f.clone(),
		}
	}

	// MulByFloat64 returns the standard scalar product of v and f.
	pub fn mul_by_f64(&self, f: f64) -> PreciseVector {
		return self.mul(prec_float(f))
	}

	// Dot returns the standard dot product of v and ov.
	pub fn dot(&self, ov: &PreciseVector) -> bigdecimal::BigDecimal{
		let a = self.clone();
		let b = ov.clone();
		return (a.x * b.x) + (a.y * b.y) + (a.z * b.z)
	}

	// Cross returns the standard cross product of v and ov.
	pub fn cross(&self, ov: PreciseVector) -> PreciseVector {
		return PreciseVector{
			x: (self.y.clone() * ov.z.clone()) - (self.z.clone() * ov.y.clone()),
			y: (self.z.clone() * ov.x.clone()) - (self.x.clone() * ov.z.clone()),
			z: (self.x.clone() * ov.y.clone()) - (self.y.clone() * ov.x.clone())
		}
	}

	// LargestComponent returns the axis that represents the largest component in this vector.
	pub fn largest_component(&self) -> r3::vector::Axis {
		let a = self.abs();
        if a.x > a.y {
            if a.x > a.z {
                r3::vector::Axis::X
            } else {
                r3::vector::Axis::Z
            }
        } else {
            if a.y > a.z {
                r3::vector::Axis::Y
            } else {
                r3::vector::Axis::Z
            }
        }
	}

	// SmallestComponent returns the axis that represents the smallest component in this vector.
	pub fn smallest_component(&self) -> r3::vector::Axis {
		let t = self.abs();
        if t.x < t.y {
            if t.x < t.z {
                r3::vector::Axis::X
            } else {
                r3::vector::Axis::Z
            }
        } else {
            if t.y < t.z {
                r3::vector::Axis::Y
            } else {
                r3::vector::Axis::Z
            }
        }
	}

}

