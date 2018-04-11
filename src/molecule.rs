// [[file:~/Workspace/Programming/gchemol/gchemol.note::942dedaa-9351-426e-9be9-cdb640ec2b75][942dedaa-9351-426e-9be9-cdb640ec2b75]]
use std::collections::HashMap;
use petgraph::prelude::*;
use petgraph as pg;

use {
    Point3D,
    Points,
    Atom,
    Bond,
};

type MolGraph = StableUnGraph<Atom, Bond>;
type AtomIndex = NodeIndex;
type BondIndex = EdgeIndex;

#[derive(Debug, Clone)]
pub struct Molecule {
    graph: MolGraph,
}

impl Molecule {
    fn new() -> Self {
        Molecule {
            graph: MolGraph::default(),
        }
    }

    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.graph.node_indices().map(move |n| &self.graph[n])
    }

    /// get positions of all atoms
    fn positions(&self) -> impl Iterator<Item = &Point3D> {
        self.atoms().map(|ref a| &a.position)
    }

    fn symbols(&self) -> impl Iterator<Item = &str> {
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

    pub fn add_bond(&mut self, atom1: AtomIndex, atom2: AtomIndex) -> BondIndex {
        let bond = Bond::default();
        self.graph.update_edge(atom1, atom2, bond)
    }
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::5052eafc-f1ab-4612-90d7-0924c3bacb16][5052eafc-f1ab-4612-90d7-0924c3bacb16]]
#[test]
fn test_molecule() {
    let mut mol = Molecule::new();
    mol.add_atom(Atom::default());
    println!("{:?}", mol);
    println!("{:?}", mol.graph);
    for a in mol.atoms() {
        println!("{:?}", a);
    }

    let a = mol.atom(0).unwrap();
    println!("{:?}", a.symbol());
}
// 5052eafc-f1ab-4612-90d7-0924c3bacb16 ends here
