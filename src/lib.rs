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
//       UPDATED:  <2018-06-04 Mon 15:37>
//===============================================================================#
// bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::53cbd3c0-e164-4bad-b535-6fd6df916650][53cbd3c0-e164-4bad-b535-6fd6df916650]]
extern crate petgraph;
#[macro_use] extern crate nom;
#[macro_use] extern crate indexmap;
extern crate nalgebra;
// extern crate cgmath;
#[macro_use] extern crate timeit;
#[macro_use] extern crate approx;
extern crate rand;
extern crate serde;
#[macro_use] extern crate serde_json;
#[macro_use] extern crate serde_derive;
extern crate itertools;
extern crate handlebars;

#[macro_use] extern crate error_chain;
// We'll put our errors in an `errors` module, and other modules in
// this crate will `use errors::*;` to get access to everything
// `error_chain!` creates.
#[macro_use]
pub mod errors {
    use std::fmt;
    use nom;

    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! {
    }

    // for nom errors
    // #[derive(Debug)]
    // pub struct NomError {
    //     desc: String,
    // }

    // impl<E: fmt::Debug + Clone> From<nom::Err<E>> for NomError {
    //     fn from(error: nom::Err<E>) -> Self {
    //         let desc = match error {
    //             nom::Err::Incomplete(needed) => format!("ran out of bytes: {:?}", needed),
    //             nom::Err::Error(context) => format!("{:?}", nom::error_to_list(&context)),
    //             nom::Err::Failure(context) => format!("{:?}", nom::error_to_list(&context)),
    //         };

    //         NomError { desc }
    //     }
    // }

    // macro_rules! nom {
    //     ($expr:expr) => {
    //         $expr.map_err(NomError::from)
    //     }
    // }
}

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
