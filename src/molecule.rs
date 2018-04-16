// [[file:~/Workspace/Programming/gchemol/gchemol.note::7e391e0e-a3e8-4c22-b881-e0425d0926bc][7e391e0e-a3e8-4c22-b881-e0425d0926bc]]
//===============================================================================#
//   DESCRIPTION:  molecular entity repsented in graph data structure
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-12 Thu 15:48>
//       UPDATED:  <2018-04-16 Mon 20:17>
//===============================================================================#
// 7e391e0e-a3e8-4c22-b881-e0425d0926bc ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::942dedaa-9351-426e-9be9-cdb640ec2b75][942dedaa-9351-426e-9be9-cdb640ec2b75]]
use std::collections::HashMap;
use petgraph;
use petgraph::prelude::*;
use std::convert;
use errors::*;

use {
    Point3D,
    Points,
    Atom,
    Bond,
};

type MolGraph = StableUnGraph<Atom, Bond>;
type AtomIndex = usize;
type BondIndex = [usize; 2];

/// Repsents any singular entity, irrespective of its nature, in order
/// to concisely express any type of chemical particle: atom,
/// molecule, ion, ion pair, radical, radical ion, complex, conformer,
/// etc.
///
/// Reference
/// ---------
/// 1. http://goldbook.iupac.org/M03986.html
/// 2. https://en.wikipedia.org/wiki/Molecular_entity
///
#[derive(Debug, Clone)]
pub struct MolecularEntity {
    /// Molecule name
    pub name: String,
    /// core data in graph
    pub graph: MolGraph,

    /// mapping atom index to NodeIndex
    atom_indices: HashMap<usize, NodeIndex>,
    /// mapping bond tuple to EdgeIndex
    bond_indices: HashMap<[usize; 2], EdgeIndex>,
}

pub type Molecule = MolecularEntity;

impl Default for Molecule {
    fn default() -> Self {
        Molecule {
            name: "default".to_string(),
            graph: MolGraph::default(),
            atom_indices: HashMap::new(),
            bond_indices: HashMap::new(),
        }
    }
}

impl Molecule {
    /// Create a new empty molecule with specific name
    pub fn new(name: &str) -> Self {
        Molecule {
            name: name.to_string(),
            ..Default::default()
        }
    }

    /// Return the number of atoms in the molecule.
    pub fn natoms(&self) -> usize {
        self.graph.node_count()
    }

    /// Return the number of bonds in the molecule.
    pub fn nbonds(&self) -> usize {
        self.graph.edge_count()
    }

    /// Construct from an existing graph
    pub fn from_graph(graph: MolGraph) -> Self{
        Molecule {
            graph: graph,
            ..Default::default()
        }
    }

    /// A convenient alias of molecular name
    pub fn title(&self) -> String {
        self.name.to_owned()
    }

    /// Return an iterator over the atoms in the molecule.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.graph.node_indices().map(move |n| &self.graph[n])
    }

    /// Return an iterator over the bonds in the molecule.
    pub fn bonds(&self) -> impl Iterator<Item = &Bond> {
        self.graph.edge_indices().map(move |e| &self.graph[e])
    }

    /// Return an iterator over positions of all atoms in the molecule.
    pub fn positions(&self) -> impl Iterator<Item = &Point3D> {
        self.atoms().map(|ref a| &a.position)
    }

    /// Return an iterator over the symbols of all  atoms in the molecule.
    pub fn symbols(&self) -> impl Iterator<Item = &str> {
        self.atoms().map(|ref a| a.symbol())
    }

    /// Set positions of atoms
    pub fn set_positions(&mut self, positions: Points)
    {
        let indices: Vec<_> = self.graph.node_indices().collect();

        for (&index, position) in indices.iter().zip(positions) {
            let mut atom = &mut self.graph[index];
            atom.position = position;
        }
    }

    /// TODO
    pub fn set_symbols(&mut self, symbols: Vec<String>) {
        //
    }
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::9924e323-dd02-49d0-ab07-41208114546f][9924e323-dd02-49d0-ab07-41208114546f]]
impl Molecule {
    /// a bunch of private methods for easy of maintaining
    /// convert atom index to graph node index
    fn get_node_index(&self, index: AtomIndex) -> Option<&NodeIndex> {
        self.atom_indices.get(&index)
    }

    /// store node index and return the associated atom index
    fn new_atom_index(&mut self, node: NodeIndex) -> AtomIndex {
        let index = self.natoms();
        // cache atom index in molecule
        self.atom_indices.insert(index, node);
        // cache atom index in atom
        let mut atom = self.get_atom_mut(index).unwrap();
        atom.index = index;

        index
    }

    /// convert a pair of atom index to graph edge index
    fn get_edge_index(&self, index1: AtomIndex, index2: AtomIndex) -> Option<&EdgeIndex> {
        self.bond_indices.get(&[index1, index2])
    }

    /// store edge index and return the assocated pair of atom indices
    fn new_bond_index(&mut self, index1: AtomIndex, index2: AtomIndex, edge: EdgeIndex) -> BondIndex {
        let mut pair = [index1, index2];
        if index1 > index2 {
            pair.swap(0, 1);
        }

        self.bond_indices.insert(pair, edge);

        pair
    }

    /// Add a single atom into molecule
    /// Return an index to atom counting from 1
    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        let i = self.graph.add_node(atom);
        self.new_atom_index(i)
    }

    /// Remove an atom from the molecule.
    /// Return the removed atom if it exists, or return None.
    pub fn remove_atom(&mut self, index: AtomIndex) -> Option<Atom> {
        if let Some(&n) = self.get_node_index(index) {
            self.graph.remove_node(n)
        } else {
            None
        }
    }

    /// TODO
    pub fn add_atoms_from(&mut self) {
        //
    }

    /// access atom by atom index
    pub fn get_atom(&self, index: AtomIndex) -> Option<&Atom> {
        if let Some(&n) = self.get_node_index(index) {
            self.graph.node_weight(n)
        } else {
            None
        }
    }

    /// mutable access to atom by atom index
    pub fn get_atom_mut(&mut self, index: AtomIndex) -> Option<&mut Atom> {
        if let Some(&n) = self.get_node_index(index) {
            self.graph.node_weight_mut(n)
        } else {
            None
        }
    }

    /// Add a single bond into molecule
    pub fn add_bond(&mut self, index1: AtomIndex, index2: AtomIndex, bond: Bond) -> Option<BondIndex> {
        if let Some(&node1) = self.get_node_index(index1) {
            if let Some(&node2) = self.get_node_index(index2) {
                // cache atom indices of the bonded pair
                let mut bond = bond;
                bond.neighbors[0] = index1;
                bond.neighbors[1] = index2;
                let e = self.graph.update_edge(node1, node2, bond);
                let b = self.new_bond_index(index1, index2, e);
                return Some(b);
            }
        };

        None
    }

    /// access bond by atom indices
    pub fn get_bond(&self, index1: AtomIndex, index2: AtomIndex) -> Option<&Bond> {
        if let Some(&e) = self.get_edge_index(index1, index2) {
            self.graph.edge_weight(e)
        } else {
            None
        }
    }

    /// Remove a bond using a pair of atom indices
    /// Return the bond if it exists, or return None
    pub fn remove_bond(&mut self, index1: AtomIndex, index2: AtomIndex) -> Option<Bond>{
        if let Some(&e) = self.get_edge_index(index1, index2) {
            self.graph.remove_edge(e)
        } else {
            None
        }
    }

    /// mutable access to bond by bond index
    pub fn get_bond_mut(&mut self, index1: AtomIndex, index2: AtomIndex) -> Option<&mut Bond> {
        if let Some(&e) = self.get_edge_index(index1, index2) {
            self.graph.edge_weight_mut(e)
        } else {
            None
        }
    }

    /// Recalculate atom indices (counting from 1)
    pub fn reorder(&mut self) {
        // reorder atom indices
        self.atom_indices.clear();
        let nodes: Vec<_> = self.graph.node_indices().collect();
        for (i, &node) in nodes.iter().enumerate() {
            let index = i + 1;
            self.atom_indices.insert(index, node);
        }

        for (i, &node) in nodes.iter().enumerate() {
            let index = i + 1;
            let atom = &mut self.graph[node];
            atom.index = index;
        }

        // reorder bond indices
        self.bond_indices.clear();
        let mut pairs = vec![];
        for e in self.graph.edge_indices() {
            let bond = &self.graph[e];
            let i = bond.neighbors[0];
            let j = bond.neighbors[1];
            if i < j {
                pairs.push((i, j, e));
            }
        }

        for (i, j, e) in pairs {
            self.new_bond_index(i, j, e);
        }
    }
}
// 9924e323-dd02-49d0-ab07-41208114546f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c936ccf3-7276-48b9-9cc7-d83b7cc257f7][c936ccf3-7276-48b9-9cc7-d83b7cc257f7]]
#[test]
fn test_molecule_remove_atom() {
    let mut mol = Molecule::default();
    let a1 = mol.add_atom(Atom::default());
    let a2 = mol.add_atom(Atom::default());
    let a3 = mol.add_atom(Atom::default());
    let a4 = mol.add_atom(Atom::default());
    let a5 = mol.add_atom(Atom::default());
    mol.get_atom_mut(a1).unwrap().position[0] = 1.0;
    mol.get_atom_mut(a2).unwrap().position[0] = 2.0;
    mol.get_atom_mut(a3).unwrap().position[0] = 2.0;
    mol.get_atom_mut(a4).unwrap().position[0] = 4.0;
    mol.get_atom_mut(a5).unwrap().position[0] = 5.0;
    // test removing atom
    {
        mol.remove_atom(a3);
        mol.add_bond(a1, a4, Bond::default());
        let a = mol.get_atom(5).unwrap();
        assert_eq!(5, a.index);
        assert_eq!(5.0, a.position[0]);
    }
    // test reorder
    {
        mol.reorder();
        let a = mol.get_atom(4).unwrap();
        assert_eq!(4, a.index);
        assert_eq!(5.0, a.position[0]);
        for b in mol.bonds() {
            println!("{:?}", b);
        }
    }
}
// c936ccf3-7276-48b9-9cc7-d83b7cc257f7 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3][a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3]]
use geometry::get_distance_matrix;
use bond::BondKind;
use data::guess_bond_kind;

impl Molecule {
    pub fn distance_matrix(&self) -> Vec<Vec<f64>>{
        let positions: Vec<_> = self.positions().map(|a| *a ).collect();
        get_distance_matrix(positions)
    }

    /// Removes all existing bonds between atoms
    pub fn unbond(&mut self) {
        self.graph.clear_edges();
    }

    /// Redefine bonds from distances based on predefined bonding lengths
    pub fn rebond(&mut self) {
        let n = self.natoms();
        let mut pairs = vec![];
        for i in 1..(n+1) {
            for j in (i+1)..(n+1) {
                let ai = self.get_atom(i).unwrap();
                let aj = self.get_atom(j).unwrap();
                let bk = guess_bond_kind(ai, aj);
                if bk != BondKind::Dummy {
                    let mut bond = Bond::default();
                    bond.kind = bk;
                    pairs.push((i, j, bond));
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
        //
    }
}
// ec7b11d2-6f13-49fd-b253-af4b213b49a3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::f0258648-03f4-41c9-949e-f3677c3b44bc][f0258648-03f4-41c9-949e-f3677c3b44bc]]
impl Molecule {
    /// fragment into a list of sub-molecules based on connectivity
    pub fn fragment(&self) -> Vec<Molecule> {
        let graph = &self.graph;

        let mut mols = vec![];
        let subgraphs = connected_component_subgraphs(graph);

        for g in subgraphs {
            mols.push(Molecule::from_graph(g));
        }

        mols
    }
}

/// Generate connected components as subgraphs
fn connected_component_subgraphs(graph: &MolGraph) -> Vec<MolGraph>{
    // get fragments from connected components
    let components = petgraph::algo::kosaraju_scc(graph);

    let mut graphs = vec![];
    for nodes in components {
        let g: MolGraph = graph.filter_map
            (
                // node closure:
                // keep atoms in the same component
                |i, n| {
                    if nodes.contains(&i) {
                        Some(n.clone())
                    } else {
                        None
                    }
                },
                // edge closure:
                // keep the bond if bonded atoms are both in the same component
                |i, e| {
                    let (n1, n2) = graph.edge_endpoints(i).unwrap();
                    if nodes.contains(&n1) && nodes.contains(&n2) {
                        Some(e.clone())
                    } else {
                        None
                    }
                }
            );
        graphs.push(g);
    }

    graphs
}
// f0258648-03f4-41c9-949e-f3677c3b44bc ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ddf54b1b-6bda-496a-8444-b9762645cc94][ddf54b1b-6bda-496a-8444-b9762645cc94]]
use std::iter::IntoIterator;
use std::fmt;

pub fn get_reduced_formula<'a, I>(symbols: I) -> String
where
    I: IntoIterator,
    I::Item: fmt::Display,
{
    let symbols: Vec<_> = symbols.into_iter()
        .map(|item| format!("{:}", item))
        .collect();

    // 1. count symbols: CCCC ==> C 4
        let mut counts = HashMap::new();
    for x in &symbols {
        let c = counts.entry(x).or_insert(0);
        *c += 1;
    }

    let mut syms: Vec<String> = Vec::new();
    let mut to_append = String::new();
    // 2. format the formula
    for (k, v) in counts {
        // 2.1 omit number if the count is 1: C1H4 ==> CH4
        let mut s = String::new();
        if v > 1 {
            s = v.to_string();
        }
        // 2.2 special treatments for C and H
        let reduced = format!("{}{}", k, s);
        if *k == "C" {
            syms.insert(0, reduced);
        } else if *k == "H" {
            to_append = reduced;
        } else {
            syms.push(reduced);
        }
    }
    // 3. final output
    syms.push(to_append);
    let formula = syms.join("");
    formula
}

impl Molecule {
    /// Return the molecule formula represented in string
    /// Return empty string if molecule containing no atom
    fn formula(&self) -> String
    {
        get_reduced_formula(self.symbols())
    }
}

#[test]
fn test_formula() {
    let symbols   = vec!["C", "H", "C", "H", "H", "H"];
    let formula = get_reduced_formula(&symbols);
    assert_eq!("C2H4", formula);
    let symbols   = vec!["C", "H", "C", "H", "H", "O", "H", "O"];
    let formula = get_reduced_formula(&symbols);
    assert_eq!("C2O2H4", formula);
}
// ddf54b1b-6bda-496a-8444-b9762645cc94 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::5052eafc-f1ab-4612-90d7-0924c3bacb16][5052eafc-f1ab-4612-90d7-0924c3bacb16]]
#[test]
fn test_molecule_basic() {
    // construct molecule
    let mut mol = Molecule::new("test");
    assert_eq!("test", mol.name);

    let mut mol = Molecule::default();
    let atom1 = Atom::new("Fe", [1.2; 3]);
    let atom2 = Atom::new("Fe", [1.0; 3]);
    let atom3 = Atom::new("C", [0.0; 3]);
    let atom4 = Atom::new("O", [2.1; 3]);
    let a1 = mol.add_atom(atom1);
    let a2 = mol.add_atom(atom2);
    let a3 = mol.add_atom(atom3);
    let a4 = mol.add_atom(atom4);
    assert_eq!(4, mol.natoms());

    let b1 = mol.add_bond(a1, a2, Bond::default());
    let b2 = mol.add_bond(a3, a4, Bond::default());
    assert_eq!(2, mol.nbonds());
    mol.remove_bond(1, 2).unwrap();
    assert_eq!(1, mol.nbonds());
    mol.add_bond(1, 4, Bond::default());
    assert_eq!(2, mol.nbonds());

    // loop over atoms
    for a in mol.atoms() {
        //
    }

    // loop over bonds
    for b in mol.bonds() {
        //
    }

    // pick a single atom
    let a = mol.get_atom(1).unwrap();
    assert_eq!("Fe", a.symbol());
    assert_eq!(1.2, a.position[0]);
    let a = mol.get_atom(a1).unwrap();
    assert_eq!("Fe", a.symbol());
    assert_eq!(1.2, a.position[0]);
}

#[test]
fn test_molecule_other() {
    let mut mol = Molecule::default();
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    //set atom positions
    let positions = [[-0.90203687,  0.62555259,  0.0081889 ],
                     [-0.54538244, -0.38325741,  0.0081889 ],
                     [-0.54536403,  1.12995078, -0.8654626 ],
                     [-1.97203687,  0.62556577,  0.0081889 ]];
    mol.set_positions(positions.to_vec());
    assert_eq!(mol.get_atom(1).unwrap().position[0], -0.90203687);

    // loop over fragments
    let frags = mol.fragment();
    for m in frags {
        m.formula();
    }
}

#[test]
fn test_molecule_rebond() {
    let atom1 = Atom::new("C", [-0.90203687,  0.62555259,  0.0081889 ]);
    let atom2 = Atom::new("H", [-0.54538244, -0.38325741,  0.0081889 ]);
    let atom3 = Atom::new("H", [-0.54536403,  1.12995078,  0.88184041]);
    let atom4 = Atom::new("H", [-0.54536403,  1.12995078, -0.8654626 ]);
    let atom5 = Atom::new("H", [-1.97203687,  0.62556577,  0.0081889 ]);

    let mut mol = Molecule::default();
    mol.add_atom(atom1);
    mol.add_atom(atom2);
    mol.add_atom(atom3);
    mol.add_atom(atom4);
    mol.add_atom(atom5);

    assert_eq!(5, mol.natoms());
    assert_eq!(0, mol.nbonds());
    mol.rebond();
    assert_eq!(4, mol.nbonds());

    use io;
    io::to_mol2file(&mol, "/tmp/test.mol2");
}
// 5052eafc-f1ab-4612-90d7-0924c3bacb16 ends here
