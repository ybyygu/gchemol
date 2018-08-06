// [[file:~/Workspace/Programming/gchemol/gchemol.note::a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3][a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3]]
use super::*;
use geometry::get_distance_matrix;
use data::guess_bond_kind;

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
}
// a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ec7b11d2-6f13-49fd-b253-af4b213b49a3][ec7b11d2-6f13-49fd-b253-af4b213b49a3]]
impl Molecule {
    /// Return an iterator of all atoms connected to a.
    fn connected() {
        unimplemented!()
    }
}
// ec7b11d2-6f13-49fd-b253-af4b213b49a3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::1f8f3a9a-dc8d-4a1e-816b-c5228ddf0607][1f8f3a9a-dc8d-4a1e-816b-c5228ddf0607]]
impl Molecule {
    // FIXME: if PBC
    pub fn distance_matrix(&self) -> Vec<Vec<f64>>{
        let positions = self.positions();
        get_distance_matrix(positions)
    }

    // TODO: improve performance
    /// Return the distance between `atom i` and `atom j`.
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
// 1f8f3a9a-dc8d-4a1e-816b-c5228ddf0607 ends here
