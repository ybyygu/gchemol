// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/geometry/geometry.note::*lib.rs][lib.rs:1]]
// for test-only uses
#[cfg(test)]
#[macro_use]
extern crate approx;

#[cfg(test)]
extern crate gchemol;

#[cfg(test)]
#[macro_use]
extern crate timeit;

mod alignment;
mod base;
mod interpolate;
mod random;
mod stats;
mod transform;

/// Re-exports important traits and types.
pub mod prelude {
    pub use super::alignment::*;
    pub use super::base::*;
    pub use super::interpolate::*;
    pub use super::transform::*;
}

/// Re-exports important functions
pub use self::random::*;

// logging, error handling, ...
mod core_utils {
    pub use guts::prelude::*;
}
// lib.rs:1 ends here
