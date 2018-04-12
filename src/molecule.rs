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

#[derive(Debug, Clone)]
pub struct Molecule {
    graph: MolGraph,
}

impl Molecule {
    /// Create a new empty molecule
    pub fn new() -> Self {
        Molecule {
            graph: MolGraph::default(),
        }
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
    fn positions(&self) -> impl Iterator<Item = &Point3D> {
        self.atoms().map(|ref a| &a.position)
    }

    /// Return an iterator over the symbols of all  atoms in the molecule.
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

    /// Add a single bond into molecule
    pub fn add_bond(&mut self, atom1: AtomIndex, atom2: AtomIndex) -> BondIndex {
        let bond = Bond::default();
        self.graph.update_edge(atom1, atom2, bond)
    }
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here

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
// 5052eafc-f1ab-4612-90d7-0924c3bacb16 ends here
