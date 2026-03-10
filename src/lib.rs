#![feature(let_chains)]
#[macro_use]
extern crate lazy_static;

#[cfg(feature = "serde")]
#[macro_use]
extern crate serde;

#[macro_use]
mod consts;

pub mod r1;
pub mod r2;
pub mod r3;

pub mod s1;

// export s2 modules directly
mod s2;
pub use crate::s2::*;
