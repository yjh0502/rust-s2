#[macro_use]
extern crate lazy_static;
extern crate float_extras;
extern crate cgmath;
extern crate rand;

#[macro_use]
mod consts;

pub mod r1;
pub mod r2;
pub mod r3;

pub mod s1;

// export s2 modules directly
mod s2;
pub use s2::*;
