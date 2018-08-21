// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::7b453b81-a6a7-41b2-b3d7-51cecd5babac][7b453b81-a6a7-41b2-b3d7-51cecd5babac]]
#![feature(test)]
extern crate test;

// for test-only uses
#[cfg(test)]
#[macro_use] extern crate approx;
#[macro_use] extern crate quicli;

extern crate nalgebra;
extern crate rand;

#[cfg(test)]
extern crate gchemol;

mod base;
mod transform;
mod random;
mod alignment;
mod interpolate;

pub mod prelude {
    pub use super::base::*;
    pub use super::transform::*;
    pub use super::random::*;
    pub use super::alignment::*;
    pub use super::interpolate::*;
}
// 7b453b81-a6a7-41b2-b3d7-51cecd5babac ends here
