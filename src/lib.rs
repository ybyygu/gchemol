// [[file:~/Workspace/Programming/gchemol/gchemol.note::bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d][bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d]]
//===============================================================================#
//   DESCRIPTION:  gchemol: a Graph-based CHEMical Objects Library
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  rewrite from scratch using rust
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 2 or upper
//       CREATED:  <2018-04-10 Tue 15:46>
//       UPDATED:  <2018-05-15 Tue 19:31>
//===============================================================================#
// bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::53cbd3c0-e164-4bad-b535-6fd6df916650][53cbd3c0-e164-4bad-b535-6fd6df916650]]
#![allow(dead_code)]

#[macro_use]
extern crate error_chain;
#[macro_use]
extern crate nom;

// We'll put our errors in an `errors` module, and other modules in
// this crate will `use errors::*;` to get access to everything
// `error_chain!` creates.
pub mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! { }
}

extern crate petgraph;
#[macro_use]
extern crate indexmap;
extern crate nalgebra;
extern crate cgmath;
#[macro_use]
extern crate timeit;
#[macro_use]
extern crate approx;
extern crate rand;
extern crate serde;
extern crate serde_json;
extern crate itertools;

pub type Point3D = [f64; 3];
pub type Points = Vec<Point3D>;

pub mod geometry;
pub mod molecule;
pub use molecule::{
    AtomKind,
    AtomKind::Element,
    AtomKind::Dummy,
    Atom,
    Bond,
    BondKind,
    Molecule};

pub mod topology;
pub mod io;
pub use io::{write_as_xyz};
pub mod data;
#[macro_use]
pub mod parser;
pub mod lattice;
pub mod formats;
// 53cbd3c0-e164-4bad-b535-6fd6df916650 ends here
