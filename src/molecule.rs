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
//       UPDATED:  <2018-04-18 Wed 21:56>
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
        let graph = MolGraph::default();
        Molecule {
            name: "default".to_string(),
            graph: graph,
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

// [[file:~/Workspace/Programming/gchemol/gchemol.note::e1d0c51a-0dd7-4977-ae54-7928ee46d373][e1d0c51a-0dd7-4977-ae54-7928ee46d373]]
use std::ops::Index;

#[derive(Debug, Clone)]
pub struct AtomsView<'a>(&'a MolGraph);

impl<'a> Index<usize> for AtomsView<'a>
{
    type Output = Atom;

    fn index(&self, index: usize) -> &Atom {
        let index = NodeIndex::new(index);
        &self.0[index]
    }
}

#[test]
fn test_atoms_view() {
    let mut mol = Molecule::default();
    mol.add_atom(Atom::new("Fe", [0.0; 3]));

    let av = AtomsView(&mol.graph);
    assert_eq!("Fe", av[0].symbol())
}
// e1d0c51a-0dd7-4977-ae54-7928ee46d373 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::5916eec2-ec7e-4525-bc6c-fade1d250a16][5916eec2-ec7e-4525-bc6c-fade1d250a16]]
#[derive(Debug, Clone)]
pub struct BondsView<'a>(&'a MolGraph);

impl<'a> Index<(usize, usize)> for BondsView<'a>
{
    type Output = Bond;

    fn index(&self, bond_index: (usize, usize)) -> &Bond {
        let i = NodeIndex::new(bond_index.0);
        let j = NodeIndex::new(bond_index.1);
        let e = self.0.find_edge(i, j).unwrap();
        &self.0[e]
    }
}

#[test]
fn test_bonds_view() {
    let mut mol = Molecule::default();
    let a1 = mol.add_atom(Atom::new("C", [0.0; 3]));
    let a2 = mol.add_atom(Atom::new("H", [1.0; 3]));
    let a3 = mol.add_atom(Atom::new("H", [2.0; 3]));
    mol.add_bond(a1, a2, Bond::default());
    let bv = BondsView(&mol.graph);
    let b = &bv[(0, 1)];
}
// 5916eec2-ec7e-4525-bc6c-fade1d250a16 ends here

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
    pub fn add_bond<T: IntoAtomIndex>(&mut self, index1: T, index2: T, bond: Bond) -> BondIndex {
        let n1 = index1.into_atom_index();
        let n2 = index2.into_atom_index();
        let e = self.graph.update_edge(n1, n2, bond);

        // cache atom indices of the bonded pair
        let mut bond = &mut self.graph[e];
        bond.index = e;

        e
    }

    /// access bond by bond index
    pub fn get_bond<T: IntoBondIndex>(&self, e: T) -> Option<&Bond> {
        let e = e.into_bond_index();
        self.graph.edge_weight(e)
    }

    /// Get the bond index between two atoms
    /// Return None if not found
    fn bond_index_between(&mut self, n1: AtomIndex, n2: AtomIndex) -> Option<BondIndex> {
        self.graph.find_edge(n1, n2)
    }

    /// Return any bond bween two atoms
    /// Return None if it does not exist
    pub fn get_bond_between<T: IntoAtomIndex>(&mut self, index1: T, index2: T) -> Option<&Bond> {
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

    /// mutable access to bond by bond index
    pub fn get_bond_mut<T: IntoBondIndex>(&mut self, b: T) -> Option<&mut Bond> {
        let b = b.into_bond_index();
        self.graph.edge_weight_mut(b)
    }

    /// Update atom indices
    pub fn reorder(&mut self) {
        let ns: Vec<_> = self.graph.node_indices().collect();
        for (i, &n) in ns.iter().enumerate() {
            let atom = &mut self.graph[n];
            atom.index = n;
        }
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

impl IntoBondIndex for Bond {
    fn into_bond_index(&self) -> BondIndex {
        self.index
    }
}
// be29e151-18c6-43cb-9586-aba0e708d38c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::72dd0c31-26e5-430b-9f67-1c5bd5220a84][72dd0c31-26e5-430b-9f67-1c5bd5220a84]]
impl Molecule {
    pub fn add_atoms_from(&mut self) {
        //
    }

    pub fn add_bonds_from(&mut self) {
        //
    }

    pub fn remove_atoms_from(&mut self) {
        //
    }

    pub fn remove_bonds_from(&mut self) {
        //
    }
}
// 72dd0c31-26e5-430b-9f67-1c5bd5220a84 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3][a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3]]
use geometry::get_distance_matrix;
use bond::BondKind;
use data::guess_bond_kind;

impl Molecule {
    pub fn distance_matrix(&self) -> Vec<Vec<f64>>{
        let positions: Vec<_> = self.positions().map(|a| *a ).collect();
        get_distance_matrix(positions)
    }

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
    pub fn unbond(&mut self) {
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
        //
    }
}
// ec7b11d2-6f13-49fd-b253-af4b213b49a3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::2a27ca30-0a99-4d5d-b544-5f5900304bbb][2a27ca30-0a99-4d5d-b544-5f5900304bbb]]
use petgraph::algo;
use geometry::euclidean_distance;
use rand::{thread_rng, Rng};

const EPSILON: f64 = 1.0E-6;

impl Molecule {
    pub fn set_position(&mut self, index: AtomIndex, position: Point3D) {
        let atom = &mut self.graph[index];
        atom.position = position;
    }

    /// Return the shortest distance numbered in bonds between two atoms
    /// Return None if them are not connected
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

    /// Translate molecule to a new location
    pub fn translate(&mut self, loc: Point3D) {
        let nodes: Vec<_> = self.graph.node_indices().collect();
        for n in nodes {
            let mut atom = &mut self.graph[n];
            for v in 0..3 {
                atom.position[v] += loc[v];
            }
        }
    }

    // return distance bounds between atoms
    // upper-tri for upper bounds
    // lower-tri for lower bounds
    fn get_distance_bounds(&self) -> Vec<Vec<f64>>{
        let mut dm = self.distance_matrix();
        // max distance between two atoms
        let max_rij = 90.0;

        let node_indices: Vec<_> = self.graph.node_indices().collect();
        let nnodes = node_indices.len();
        for i in 0..nnodes {
            for j in (i+1)..nnodes {
                let node_i = node_indices[i];
                let node_j = node_indices[j];
                let atom_i = &self.graph[node_i];
                let atom_j = &self.graph[node_j];

                // use vdw radii as the lower bound for non-bonded pair
                let vri = atom_i.vdw_radius().unwrap();
                let vrj = atom_j.vdw_radius().unwrap();
                let vrij = vri + vrj;

                // use covalent radii as the lower bound for bonded pair
                let cri = atom_i.covalent_radius().unwrap();
                let crj = atom_j.covalent_radius().unwrap();
                let crij = cri + crj;

                // make sure vdw radius larger than covalent radius (usually it is)
                let mut bound = [crij, vrij];
                if crij > vrij {
                    bound.swap(0, 1);
                }

                let dij = dm[i][j];
                // if i and j is directly bonded
                // set covalent radius as the lower bound
                // or set vdw radius as the lower bound if not bonded
                if let Some(nb) = self.nbonds_between(node_i, node_j) {
                    if nb == 1 {
                        dm[i][j] = bound[0];
                        dm[j][i] = bound[0];
                    } else if nb == 2 {
                        if dij > bound[1] && dij < max_rij {
                            dm[i][j] = bound[1];
                            dm[j][i] = max_rij;
                        } else {
                            dm[i][j] = bound[1];
                            dm[j][i] = max_rij;
                        }
                    } else {
                        if dij > bound[1] && dij < max_rij {
                            dm[i][j] = bound[1];
                            dm[j][i] = max_rij;
                        } else {
                            dm[i][j] = bound[1];
                            dm[j][i] = max_rij;
                        }
                    }
                    println!("{:?}", ((node_i.index(), node_j.index()), bound, dij, nb, dm[i][j], dm[j][i]));
                } else {
                    dm[i][j] = bound[1];
                    dm[j][i] = max_rij;
                }
            }
        }

        dm
    }

    /// Clean up molecule geometry
    pub fn clean(&mut self) {
        let mut rng = thread_rng();

        // parameters for SPE
        let maxlam = 1.0;
        let minlam = 0.001;
        let maxcycle = 500;
        let maxstep = 100;
        let dlam = (maxlam - minlam)/(maxcycle as f64);
        let mut lam = maxlam;

        let indices: Vec<_> = self.graph.node_indices().collect();
        let orign_distances = self.distance_matrix();
        let bounds = self.get_distance_bounds();

        // enter spe loop
        for icycle in 0..maxcycle {
            lam -= dlam;
            // choose the first random atom
            let &node_i = rng.choose(&indices).expect("node_i");
            for istep in 0..maxstep {
                println!("cycle {:}, step {:}", icycle, istep);
                // choose the second random atom
                let &node_j = rng.choose(&indices).expect("node_j");
                let lij = bounds[node_i.index()][node_j.index()];
                let uij = bounds[node_j.index()][node_i.index()];
                let mut bound = [lij, uij];
                if lij > uij {
                    bound.swap(0, 1);
                }
                let lij = bound[0];
                let uij = bound[1];

                // calculate current distance between node_i and node_j
                let pi = self.get_atom(node_i).expect("atom i from node_i").position;
                let pj = self.get_atom(node_j).expect("atom j from node_j").position;
                let dij = euclidean_distance(pi, pj);

                // update coordinates or not
                if dij >= lij && dij <= uij {
                    continue;
                }

                let oij = orign_distances[node_i.index()][node_j.index()];
                let mut rij = lij;
                if oij >= lij && oij <= uij {
                    rij = oij;
                }

                {
                    let mut bb = [node_i.index(), node_j.index()];
                    if bb[0] > bb[1] {
                        bb.swap(0, 1);
                    }
                    println!("{:?}", (bb, dij, rij, bound));
                }

                // calculate position shift
                let mut pi = pi;
                let mut pj = pj;
                let disparity = lam*(rij - dij) / (dij + EPSILON);
                for v in 0..3 {
                    pi[v] += 0.5 * disparity * (pi[v] - pj[v]);
                    pj[v] += 0.5 * disparity * (pj[v] - pi[v]);
                }

                self.set_position(node_i, pi);
                self.set_position(node_j, pj);
            }
        }
        println!("{:?}", bounds);
    }
}

#[test]
fn test_molecule_clean() {
    let mut mol = Molecule::from_file("/tmp/test.mol2").unwrap();
    mol.clean();
    mol.to_file("/tmp/test2.mol2");
}
// 2a27ca30-0a99-4d5d-b544-5f5900304bbb ends here

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
    mol.remove_bond_between(a1, a2).expect("failed to remove bond between a1 and a2");
    mol.remove_bond(b2).expect("failed to remove bond b2");
    assert_eq!(0, mol.nbonds());
    let b14 = mol.add_bond(a1, a4, Bond::default());
    assert_eq!(1, mol.nbonds());

    // bonded partners
    let real_b14 = mol.get_bond(b14).expect("failed to get bond b14");
    assert_eq!(real_b14.index, b14);
    let (n1, n4) = mol.partners(&b14).expect("failed to get bond partners using bond index");
    assert_eq!(n1.index(), a1.index());
    assert_eq!(n4.index(), a4.index());
    // get partners using bond
    let (n1, n4) = mol.partners(real_b14).expect("failed to get bond partners using bond struct");
    assert_eq!(n1.index(), a1.index());
    assert_eq!(n4.index(), a4.index());

    // get atom neighbors
    let indices = mol.neighbors(a1);
    assert_eq!(1, indices.len());
    assert!(indices.contains(&a4));

    // loop over atoms
    for a in mol.atoms() {
        //
    }

    // loop over bonds
    for b in mol.bonds() {
        //
    }

    // pick a single atom
    let a = mol.get_atom(0).expect("failed to get atom with index 0");
    assert_eq!("Fe", a.symbol());
    assert_eq!(1.2, a.position[0]);
    let a = mol.get_atom(a1).expect("failed to get atom a1");
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
    let a = mol.get_atom(0).expect("failed to get atom with index 0");
    assert_eq!(a.position[0], -0.90203687);

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
}
// 5052eafc-f1ab-4612-90d7-0924c3bacb16 ends here
