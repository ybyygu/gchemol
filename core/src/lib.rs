// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*lib.rs][lib.rs:1]]
extern crate petgraph;
#[macro_use] extern crate indexmap;
extern crate nalgebra;
extern crate serde;
#[macro_use] extern crate serde_json;
#[macro_use] extern crate serde_derive;
extern crate itertools;

// for tests
#[cfg(test)]
#[macro_use] extern crate approx;

pub type Point3D = [f64; 3];
pub type Points = Vec<Point3D>;

extern crate gchemol_geometry;
pub mod geometry {
    pub use gchemol_geometry::*;
}

pub mod molecule;
pub mod lattice;
pub mod data;

pub use crate::molecule::{
    Atom,
    AtomKind,
    Bond,
    BondKind,
    Molecule,
    AtomIndex,
    BondIndex,
};

pub use crate::lattice::{
    Lattice,
    Supercell,
};

mod core_utils {
    pub use quicli::prelude::*;
    pub type Result<T> = ::std::result::Result<T, Error>;
}
// lib.rs:1 ends here
