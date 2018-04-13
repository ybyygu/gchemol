// [[file:~/Workspace/Programming/gchemol/gchemol.note::942dedaa-9351-426e-9be9-cdb640ec2b75][942dedaa-9351-426e-9be9-cdb640ec2b75]]
use std::collections::HashMap;
use petgraph::prelude::*;
use petgraph as pg;
use errors::*;

use {
    Point3D,
    Points,
    Atom,
    Bond,
};

type MolGraph = StableUnGraph<Atom, Bond>;
type AtomIndex = NodeIndex;
type BondIndex = EdgeIndex;

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
    graph: MolGraph,
}

pub type Molecule = MolecularEntity;

impl Molecule {
    /// Create a new empty molecule
    pub fn new() -> Self {
        Molecule {
            graph: MolGraph::default(),
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
            graph: graph
        }
    }

    /// Return an iterator over the atoms in the molecule.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.graph.node_indices().map(move |n| &self.graph[n])
    }

    /// Return an iterator over positions of all atoms in the molecule.
    pub fn positions(&self) -> impl Iterator<Item = &Point3D> {
        self.atoms().map(|ref a| &a.position)
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

    /// Return an iterator over the symbols of all  atoms in the molecule.
    pub fn symbols(&self) -> impl Iterator<Item = &str> {
        self.atoms().map(|ref a| a.symbol())
    }

    /// access atom by atom index
    pub fn atom(&self, index: usize) -> Option<&Atom> {
        let i = AtomIndex::new(index);

        self.graph.node_weight(i)
    }

    /// Add a single atom into molecule
    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        self.graph.add_node(atom)
    }

    /// Remove an atom from the molecule.
    /// Return atom if it exists, or return None.
    pub fn remove_atom(&mut self, a: AtomIndex) -> Option<Atom> {
        self.graph.remove_node(a)
    }
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::be29e151-18c6-43cb-9586-aba0e708d38c][be29e151-18c6-43cb-9586-aba0e708d38c]]
pub trait IntoAtomIndex {
    fn into_atom_index(&self) -> AtomIndex;
}

impl IntoAtomIndex for usize {
    fn into_atom_index(&self) -> AtomIndex {
        // atom index counting from zero
        debug_assert!(*self > 0);
        AtomIndex::new(*self - 1)
    }
}

impl IntoAtomIndex for AtomIndex {
    fn into_atom_index(&self) -> AtomIndex {
        *self
    }
}
// be29e151-18c6-43cb-9586-aba0e708d38c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3][a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3]]
use geometry::get_distance_matrix;

impl Molecule {
    fn distance_matrix(&self) -> Vec<Vec<f64>>{
        let positions: Vec<_> = self.positions().map(|a| *a ).collect();
        get_distance_matrix(positions)
    }
}
// a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::66630db1-08e3-479e-b59f-00c5c3b08164][66630db1-08e3-479e-b59f-00c5c3b08164]]
impl Molecule {
    /// Add a single bond into molecule
    pub fn add_bond(&mut self, atom1: AtomIndex, atom2: AtomIndex) -> BondIndex {
        let bond = Bond::default();
        self.graph.update_edge(atom1, atom2, bond)
    }

    /// Remove a bond using atom indices
    /// Return the bond if it exists, or return None
    pub fn remove_bond_between<T: IntoAtomIndex>(&mut self, a: T, b: T) -> Option<Bond> {
        let a = a.into_atom_index();
        let b = b.into_atom_index();
        if let Some(e) = self.graph.find_edge(a, b) {
            self.remove_bond(e)
        } else {
            None
        }
    }

    /// Remove a bond using bond index
    /// Return the bond if it exists, or return None
    pub fn remove_bond(&mut self, e: BondIndex) -> Option<Bond>{
        self.graph.remove_edge(e)
    }

}
// 66630db1-08e3-479e-b59f-00c5c3b08164 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::f556412b-874d-4c2c-ba79-9fb3edbefae1][f556412b-874d-4c2c-ba79-9fb3edbefae1]]
use bond::BondKind;
use data::guess_bond_kind;

impl Molecule {
    /// Removes all existing bonds between atoms
    pub fn unbond(&mut self) {
        self.graph.clear_edges();
    }

    /// Redefine bonds from distances based on predefined bonding lengths
    pub fn rebond(&mut self) {
        let indices: Vec<_> = self.graph.node_indices().collect();
        let n = indices.len();
        let mut pairs = vec![];
        for i in 0..n {
            for j in (i+1)..n {
                let vi = indices[i];
                let vj = indices[j];
                let ai = &self.graph[vi];
                let aj = &self.graph[vj];
                let bk = guess_bond_kind(ai, aj);
                if bk != BondKind::Dummy {
                    pairs.push((vi, vj));
                }
            }
        }

        for (vi, vj) in pairs {
            self.add_bond(vi, vj);
        }
    }
}
// f556412b-874d-4c2c-ba79-9fb3edbefae1 ends here

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
    fn fragment(&self) -> Vec<Molecule> {
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
    let components = pg::algo::kosaraju_scc(graph);

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
fn test_molecule() {
    let mut mol = Molecule::new();
    let atom1 = Atom::new("Fe", [1.2; 3]);
    let atom2 = Atom::new("Fe", [1.0; 3]);
    let atom3 = Atom::new("C", [0.0; 3]);
    let atom4 = Atom::new("O", [2.1; 3]);
    let a1 = mol.add_atom(atom1);
    let a2 = mol.add_atom(atom2);
    let a3 = mol.add_atom(atom3);
    let a4 = mol.add_atom(atom4);

    let b1 = mol.add_bond(a1, a2);
    let b2 = mol.add_bond(a3, a4);

    // loop over atoms
    for a in mol.atoms() {
        println!("got {:?}", a);
    }

    //set atom positions
    let positions = [[-0.90203687,  0.62555259,  0.0081889 ],
                     [-0.54538244, -0.38325741,  0.0081889 ],
                     [-0.54536403,  1.12995078, -0.8654626 ],
                     [-1.97203687,  0.62556577,  0.0081889 ]];
    mol.set_positions(positions.to_vec());
    assert_eq!(mol.atom(0).unwrap().position[0], -0.90203687);

    // pick a single atom
    let a = mol.atom(0).unwrap();
    println!("{:?}", a.symbol());
    println!("{:?}", mol.formula());

    // loop over fragments
    let frags = mol.fragment();
    for m in frags {
        println!("{:?}", m.formula());
    }
}

#[test]
fn test_molecule_rebond() {
    let atom1 = Atom::new("C", [-0.90203687,  0.62555259,  0.0081889 ]);
    let atom2 = Atom::new("H", [-0.54538244, -0.38325741,  0.0081889 ]);
    let atom3 = Atom::new("H", [-0.54536403,  1.12995078,  0.88184041]);
    let atom4 = Atom::new("H", [-0.54536403,  1.12995078, -0.8654626 ]);
    let atom5 = Atom::new("H", [-1.97203687,  0.62556577,  0.0081889 ]);

    let mut mol = Molecule::new();
    mol.add_atom(atom1);
    mol.add_atom(atom2);
    mol.add_atom(atom3);
    mol.add_atom(atom4);
    mol.add_atom(atom5);

    assert_eq!(5, mol.natoms());
    assert_eq!(0, mol.nbonds());
    mol.rebond();
    assert_eq!(4, mol.nbonds());

    let x = mol.remove_bond_between(1, 4);
    assert_eq!(3, mol.nbonds());
}
// 5052eafc-f1ab-4612-90d7-0924c3bacb16 ends here
