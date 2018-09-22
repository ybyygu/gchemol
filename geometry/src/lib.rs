// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::*lib.rs][lib.rs:1]]
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

#[cfg(test)]
#[macro_use] extern crate timeit;

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
// lib.rs:1 ends here
