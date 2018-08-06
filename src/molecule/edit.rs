// [[file:~/Workspace/Programming/gchemol/gchemol.note::548b6f53-eb69-42e5-bc84-d4e52cc17888][548b6f53-eb69-42e5-bc84-d4e52cc17888]]
use super::*;
// 548b6f53-eb69-42e5-bc84-d4e52cc17888 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::9924e323-dd02-49d0-ab07-41208114546f][9924e323-dd02-49d0-ab07-41208114546f]]
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
    pub fn remove_atom<T: IntoAtomIndex>(&mut self, index: T) -> Option<Atom> {
        let n = index.into_atom_index();
        self.graph.remove_node(n)
    }

    /// access atom by atom index
    pub fn get_atom<T: IntoAtomIndex>(&self, index: T) -> Option<&Atom> {
        let n = index.into_atom_index();
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
    pub fn remove_bond<T: IntoBondIndex>(&mut self, b: T) -> Option<Bond>{
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
// 9924e323-dd02-49d0-ab07-41208114546f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::be29e151-18c6-43cb-9586-aba0e708d38c][be29e151-18c6-43cb-9586-aba0e708d38c]]
pub trait IntoAtomIndex {
    fn into_atom_index(&self) -> AtomIndex;
}

impl IntoAtomIndex for usize {
    fn into_atom_index(&self) -> AtomIndex {
        // atom index counting from zero
        AtomIndex::new(*self)
    }
}

impl IntoAtomIndex for AtomIndex {
    fn into_atom_index(&self) -> AtomIndex {
        *self
    }
}

pub trait IntoBondIndex {
    fn into_bond_index(&self) -> BondIndex;
}

impl IntoBondIndex for usize {
    fn into_bond_index(&self) -> BondIndex {
        // atom index counting from zero
        BondIndex::new(*self)
    }
}

impl IntoBondIndex for BondIndex {
    fn into_bond_index(&self) -> BondIndex {
        *self
    }
}
// be29e151-18c6-43cb-9586-aba0e708d38c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::72dd0c31-26e5-430b-9f67-1c5bd5220a84][72dd0c31-26e5-430b-9f67-1c5bd5220a84]]
impl Molecule {
    /// add many atoms from a hashmap
    pub fn add_atoms_from(&mut self, atoms: HashMap<&str, Atom>) -> Result<()>{
        unimplemented!()
    }

    pub fn add_bonds_from(&mut self, bonds: HashMap<(String, String), Bond>) -> Result<()>{
        unimplemented!()
    }

    pub fn remove_atoms_from(&mut self) {
        unimplemented!()
    }

    pub fn remove_bonds_from(&mut self) {
        unimplemented!()
    }
}
// 72dd0c31-26e5-430b-9f67-1c5bd5220a84 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::2a27ca30-0a99-4d5d-b544-5f5900304bbb][2a27ca30-0a99-4d5d-b544-5f5900304bbb]]
use petgraph::algo;
use rand::{thread_rng, Rng};

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

#[test]
fn test_molecule_center() {
    let mol = Molecule::from_file("tests/files/mol2/alanine-gv.mol2").expect("mol2 gv");
    let pc = mol.center_of_geometry();
    let pe = [-2.31413333, -1.24455833,  0.41005833];

    for i in 0..3 {
        assert_relative_eq!(pe[i], pc[i], epsilon=1e-4);
    }

    // TODO
    // center of mass
    // let pe = [-1.8295483 , -1.10700382,  0.25332597];
}
// 2a27ca30-0a99-4d5d-b544-5f5900304bbb ends here
