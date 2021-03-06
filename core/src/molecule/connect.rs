// bonding
// chemical bonds related methods

use super::*;
use crate::data::guess_bond_kind;

impl Molecule {
    /// Access the bonded atom indices for bond `b`
    pub fn partners<T: IntoBondIndex>(&self, b: &T) -> Option<(AtomIndex, AtomIndex)>{
        let b = b.into_bond_index();
        self.graph.edge_endpoints(b)
    }

    /// Return all connected atoms with `a`
    pub fn neighbors(&self, a: AtomIndex) -> Vec<AtomIndex> {
        self.graph.neighbors(a).collect()
    }

    /// Removes all existing bonds between atoms
    pub fn unbound(&mut self) {
        self.graph.clear_edges();
    }

    /// Redefine bonds from distances based on predefined bonding lengths
    pub fn rebond(&mut self) {
        if self.lattice.is_some() {
            unimplemented!();
        }

        let n = self.natoms();
        let mut pairs = vec![];
        for ni in self.graph.node_indices() {
            for nj in self.graph.node_indices() {
                if ni < nj {
                    let ai = &self.graph[ni];
                    let aj = &self.graph[nj];
                    let bk = guess_bond_kind(ai, aj);
                    if bk != BondKind::Dummy {
                        let mut bond = Bond::default();
                        bond.kind = bk;
                        pairs.push((ni, nj, bond));
                    }
                }
            }
        }

        for (i, j, b) in pairs {
            self.add_bond(i, j, b);
        }
    }

    /// Return an iterator of all atoms connected to a.
    fn connected(&self) {
        unimplemented!()
    }
}

// distance

use crate::geometry::prelude::*;

impl Molecule {
    /// Return the distance matrix (nalgebra dmatrix). Will apply mic if
    /// possible.
    pub fn distance_matrix(&self) -> DMatrixf {
        let positions = self.positions();
        let mut dm = positions.distance_matrix();

        // apply mic for periodic structure
        let natoms = positions.len();
        if let Some(mut lat) = self.lattice {
            for i in 0..natoms {
                for j in (i+1)..natoms {
                    let pi = positions[i];
                    let pj = positions[j];
                    let dij = lat.distance(pi, pj);
                    dm[(i, j)] = dij;
                    dm[(j, i)] = dij;
                }
            }
        }

        dm
    }

    /// Return the distance between `atom i` and `atom j`.
    /// Return None if atom index out of range.
    ///
    /// Force periodic structure, this method will return the distance under the
    /// minimum image convention.
    pub fn distance(&self, i: AtomIndex, j: AtomIndex) -> Option<f64> {
        if let Some(ai) = self.get_atom(i) {
            if let Some(aj) = self.get_atom(j) {
                if let Some(mut lat) = self.lattice {
                    let pi = ai.position();
                    let pj = aj.position();
                    let dij = lat.distance(pi, pj);
                    Some(dij)
                } else {
                    Some(ai.distance(aj))
                }
            } else {
                None
            }
        } else {
            None
        }
    }
}
