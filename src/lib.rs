// [[file:~/Workspace/Programming/gchemol/gchemol.note::53cbd3c0-e164-4bad-b535-6fd6df916650][53cbd3c0-e164-4bad-b535-6fd6df916650]]
#![feature(conservative_impl_trait)]
#![allow(dead_code)]

#[macro_use]
extern crate error_chain;

// We'll put our errors in an `errors` module, and other modules in
// this crate will `use errors::*;` to get access to everything
// `error_chain!` creates.
pub mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! { }
}

extern crate petgraph;
extern crate cgmath;
#[macro_use]
extern crate timeit;
#[macro_use]
extern crate approx;
extern crate rand;
extern crate itertools;

pub type Point3D = [f64; 3];
pub type Points = Vec<Point3D>;

pub mod geometry;
pub mod element;
pub use element::AtomKind;
pub use element::AtomKind::{Element, Dummy};
pub mod atom;
pub use atom::Atom;
pub mod bond;
pub use bond::Bond;
pub mod molecule;
pub mod topology;
pub mod io;
pub use io::write_as_xyz;
pub mod data;
// 53cbd3c0-e164-4bad-b535-6fd6df916650 ends here
