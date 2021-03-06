// base

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*base][base:1]]
use super::*;

impl Molecule {
    // FIXME: opt performance
    /// Set positions of atoms
    pub fn set_positions(&mut self, positions: &[[f64; 3]]) -> Result<()>
    {
        let indices = self.sites();
        if indices.len() != positions.len() {
            bail!("the number of cartesian coordinates is different from the number of atoms in molecule.")
        }

        for (&index, &position) in indices.iter().zip(positions) {
            let mut atom = &mut self.graph[index];
            atom.set_position(position);
        }

        Ok(())
    }

    /// Set element symbols
    pub fn set_symbols(&mut self, symbols: &[&str]) -> Result<()> {
        let indices = self.sites();

        for (&index, symbol) in indices.iter().zip(symbols) {
            let mut atom = &mut self.graph[index];
            atom.set_symbol(symbol);
        }

        Ok(())
    }
}
// base:1 ends here

// one by one

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*one by one][one by one:1]]
impl Molecule {
    /// Add a single atom into molecule
    /// Return an index to atom counting from 1
    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        let n = self.graph.add_node(atom);
        let atom = self.get_atom_mut(n).unwrap();
        atom.index = n;

        n
    }

    /// Remove an atom from the molecule.
    /// Return the removed atom if it exists, or return None.
    pub fn remove_atom(&mut self, n: AtomIndex) -> Option<Atom> {
        self.graph.remove_node(n)
    }

    /// access atom by atom index
    pub fn get_atom(&self, n: AtomIndex) -> Option<&Atom> {
        self.graph.node_weight(n)
    }

    /// mutable access to atom by atom index
    pub fn get_atom_mut<T: IntoAtomIndex>(&mut self, index: T) -> Option<&mut Atom> {
        let n = index.into_atom_index();
        self.graph.node_weight_mut(n)
    }

    /// Add a single bond into molecule specified by atom indices
    /// Will panic if corresponding atoms does not exist
    /// The existing bond data will be replaced if n1 already bonded with n2
    /// Return a bond index pointing to bond data
    pub fn add_bond(&mut self, n1: AtomIndex, n2: AtomIndex, bond: Bond) -> BondIndex {
        let e = self.graph.update_edge(n1, n2, bond);

        // cache the pair of atoms
        let mut bond = &mut self.graph[e];
        bond.index = e;

        e
    }

    /// Access bond by bond index
    pub fn get_bond<T: IntoBondIndex>(&self, e: T) -> Option<&Bond> {
        let e = e.into_bond_index();
        self.graph.edge_weight(e)
    }

    /// Get the bond index between two atoms.
    /// Return None if not found
    fn bond_index_between(&self, n1: AtomIndex, n2: AtomIndex) -> Option<BondIndex> {
        self.graph.find_edge(n1, n2)
    }

    /// Return any bond bween two atoms.
    /// Return None if it does not exist
    pub fn get_bond_between<T: IntoAtomIndex>(&self, index1: T, index2: T) -> Option<&Bond> {
        let n1 = index1.into_atom_index();
        let n2 = index2.into_atom_index();

        if let Some(e) = self.bond_index_between(n1, n2) {
            self.get_bond(e)
        } else {
            None
        }
    }

    /// Remove a bond specified by bond index
    /// Return the removed bond if it exists, or return None
    pub fn remove_bond<T: IntoBondIndex>(&mut self, b: T) -> Option<Bond> {
        let b = b.into_bond_index();
        self.graph.remove_edge(b)
    }

    /// Remove any bond bween two atoms
    /// Return the removed bond if it exists, or return None
    pub fn remove_bond_between<T: IntoAtomIndex>(&mut self, index1: T, index2: T) -> Option<Bond> {
        let n1 = index1.into_atom_index();
        let n2 = index2.into_atom_index();

        if let Some(e) = self.bond_index_between(n1, n2) {
            self.remove_bond(e)
        } else {
            None
        }
    }

    /// Removes all bonds between two selections to respect pymol's unbond command.
    ///
    /// Parameters
    /// ----------
    /// atom_indices1: the first collection of atoms
    ///
    /// atom_indices2: the other collection of atoms
    ///
    /// Reference
    /// ---------
    /// https://pymolwiki.org/index.php/Unbond
    ///
    pub fn unbond(&mut self, atom_indices1: Vec<AtomIndex>, atom_indices2: Vec<AtomIndex>) {
        for &index1 in atom_indices1.iter() {
            for &index2 in atom_indices2.iter() {
                self.remove_bond_between(index1, index2);
            }
        }
    }

    /// mutable access to bond by bond index
    pub fn get_bond_mut<T: IntoBondIndex>(&mut self, b: T) -> Option<&mut Bond> {
        let b = b.into_bond_index();
        self.graph.edge_weight_mut(b)
    }
}
// one by one:1 ends here

impl Molecule {
    /// add many atoms from a hashmap
    pub fn add_atoms_from<T>(&mut self, atoms: T)
    where
        T: IntoIterator,
        T::Item: Into<Atom>,
    {
        for a in atoms {
            self.add_atom(a.into());
        }
    }

    pub fn add_bonds_from(&mut self, bonds: HashMap<(String, String), Bond>) -> Result<()> {
        unimplemented!()
    }

    pub fn remove_atoms_from(&mut self) {
        unimplemented!()
    }

    pub fn remove_bonds_from(&mut self) {
        unimplemented!()
    }
}

// translation

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*translation][translation:1]]
use petgraph::algo;

const EPSILON: f64 = 1.0E-6;

impl Molecule {
    /// Set atom position
    pub fn set_position(&mut self, index: AtomIndex, position: Point3D) {
        let atom = &mut self.graph[index];
        atom.set_position(position);
    }

    /// Return the shortest distance numbered in bonds between two atoms
    /// Return None if they are not connected
    pub fn nbonds_between(&self, index1: AtomIndex, index2: AtomIndex) -> Option<usize> {
        let path = algo::astar(&self.graph,
                               index1,
                               |finish| finish == index2,
                               |e| 1,
                               |_| 0);
        if let Some((n, _)) = path {
            Some(n)
        } else {
            None
        }
    }

    /// Return the shortest path.
    /// Return None if them are not connected
    pub fn path_between(&self, index1: AtomIndex, index2: AtomIndex) -> Option<Vec<AtomIndex>> {
        let path = algo::astar(&self.graph,
                               index1,
                               |finish| finish == index2,
                               |e| 1,
                               |_| 0);
        if let Some((_, p)) = path {
            Some(p)
        } else {
            None
        }
    }

    /// Translate atomic positions by a displacement
    pub fn translate(&mut self, displacement: Point3D) {
        let nodes: Vec<_> = self.graph.node_indices().collect();
        for n in nodes {
            let mut atom = &mut self.graph[n];
            let mut position = atom.position();
            for v in 0..3 {
                position[v] += displacement[v];
            }
            atom.set_position(position);
        }
    }

    /// Return the center of mass of molecule (COM).
    pub fn center_of_mass(&self) -> Point3D {
        unimplemented!()
    }

    /// Return the center of geometry of molecule (COG).
    pub fn center_of_geometry(&self) -> Point3D {
        let mut p = [0.0; 3];
        for [x, y, z] in self.positions() {
            p[0] += x;
            p[1] += y;
            p[2] += z;
        }

        let n = self.natoms() as f64;
        p[0] /= n;
        p[1] /= n;
        p[2] /= n;

        p
    }

    /// Center the molecule around its center of geometry
    pub fn recenter(&mut self) {
        let mut p = self.center_of_geometry();
        for i in 0..3 {
            p[i] *= -1.0;
        }
        self.translate(p);
    }
}
// translation:1 ends here
