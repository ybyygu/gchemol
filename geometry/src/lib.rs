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

#[cfg(test)]
extern crate gchemol;

#[cfg(test)]
#[macro_use] extern crate timeit;

mod base;
mod transform;
mod random;
mod alignment;
mod interpolate;

/// Re-exports important traits and types.
pub mod prelude {
    pub use super::base::*;
    pub use super::transform::*;
    pub use super::alignment::*;
    pub use super::interpolate::*;
}

/// Re-exports important functions
pub use self::random::*;
// lib.rs:1 ends here
