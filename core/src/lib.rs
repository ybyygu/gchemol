// [[file:~/Workspace/Programming/gchemol/core/gchemol-core.note::42fb78a1-aae3-40aa-8ae0-aa68b10f3f7e][42fb78a1-aae3-40aa-8ae0-aa68b10f3f7e]]
extern crate petgraph;
#[macro_use] extern crate indexmap;
extern crate nalgebra;
extern crate serde;
#[macro_use] extern crate serde_json;
#[macro_use] extern crate serde_derive;
extern crate itertools;
#[macro_use] extern crate quicli;

// for test-only uses
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

pub use molecule::{
    Atom,
    Bond,
    BondKind,
    Molecule,
    AtomIndex,
    BondIndex,
};

pub use lattice::{
    Lattice,
};
// 42fb78a1-aae3-40aa-8ae0-aa68b10f3f7e ends here
