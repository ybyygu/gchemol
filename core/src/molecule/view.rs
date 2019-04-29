// atoms view

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*atoms%20view][atoms view:1]]
use std::ops::Index;
use super::*;

/// A list-like object providing a convenient view on atoms in molecule
#[derive(Debug, Clone)]
pub struct AtomsView<'a> {
    /// mapping a positive integer to internal graph node index
    mapping: HashMap<usize, AtomIndex>,
    /// parent molecule struct
    parent: &'a Molecule,

    // current position in iteration
    cur: usize,
}

impl<'a> AtomsView<'a> {
    pub fn new(mol: &'a Molecule) -> Self {
        // use a hash map to cache graph node indices
        let mut mapping = HashMap::new();
        let mut i = 1;
        for index in mol.graph.node_indices() {
            mapping.insert(i, index);
            i += 1;
        }

        AtomsView {
            mapping,
            parent: mol,
            cur: 0,
        }
    }
}

/// Index the atoms in `Molecule` by index counting from 1
/// Will panic if index is invalid
impl<'a> Index<usize> for AtomsView<'a>
{
    type Output = Atom;

    fn index(&self, index: usize) -> &Atom {
        let n = self.mapping[&index];
        &self.parent.graph[n]
    }
}

impl<'a> Iterator for AtomsView<'a> {
    type Item = (usize, &'a Atom);

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur >= self.mapping.len() {
            None
        } else {
            self.cur += 1;
            let n = self.mapping[&self.cur];
            let a = &self.parent.graph[n];

            Some((self.cur, &a))
        }
    }
}

#[test]
fn test_atom_view() {
    let mut mol = Molecule::default();
    mol.add_atom(Atom::new("Fe", [0.0; 3]));
    mol.add_atom(Atom::new("C", [0.0; 3]));

    let av = AtomsView::new(&mol);
    assert_eq!("Fe", av[1].symbol());

    // iterate with a index (counting from 1) and an atom object
    for (i, a) in av {
        //
    }
}

impl Molecule {
    pub fn view_atoms(&self) -> AtomsView {
        AtomsView::new(&self)
    }
}
// atoms view:1 ends here

// bonds view

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*bonds%20view][bonds view:1]]
use indexmap::{IndexMap, indexmap};

/// A list-like object providing a convenient view on bonds in molecule
#[derive(Debug, Clone)]
pub struct BondsView<'a> {
    mapping: IndexMap<(usize, usize), BondIndex>,
    parent: &'a Molecule,

    // current position in iteration
    cur: usize,
}

impl<'a> BondsView<'a> {
    pub fn new(mol: &'a Molecule) -> Self {
        // reverse mapping atom id and internal graph node index
        let mut d = HashMap::new();
        let mut i = 1;
        for n in mol.graph.node_indices() {
            d.insert(n, i);
            i += 1;
        }
        // use a hash map to cache graph edge indices
        let mut mapping = indexmap!{};
        for e in mol.graph.edge_indices() {
            let (ni, nj) = mol.graph.edge_endpoints(e).expect("bondview endpoints");
            let ai = d[&ni];
            let aj = d[&nj];
            // make sure ai is always smaller than aj
            if ai < aj {
                mapping.insert((ai, aj), e);
            } else {
                mapping.insert((aj, ai), e);
            }
        }

        BondsView {
            mapping,
            parent: mol,
            cur: 0,
        }
    }
}

impl<'a> Index<(usize, usize)> for BondsView<'a>
{
    type Output = Bond;

    fn index(&self, b: (usize, usize)) -> &Bond {
        // make sure the first index number is always smaller
        let e = if b.0 < b.1 {
            self.mapping[&b]
        } else {
            self.mapping[&(b.1, b.0)]
        };

        &self.parent.graph[e]
    }
}

impl<'a> Iterator for BondsView<'a> {
    type Item = (usize, usize, &'a Bond);

    /// return a tuple in (index_i, index_j, bond)
    /// index_i is always smaller than index_j
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur >= self.mapping.len() {
            None
        } else {
            let (&(i, j), &e) = self.mapping.get_index(self.cur).expect("bondview: get bond by index");
            let b = &self.parent.get_bond(e).expect("bondview: get bond by edge index");
            self.cur += 1;

            if i < j {
                Some((i, j, &b))
            } else {
                Some((j, i, &b))
            }
        }
    }
}

#[test]
fn test_bonds_view() {
    let mut mol = Molecule::default();
    let a1 = mol.add_atom(Atom::new("C", [0.0; 3]));
    let a2 = mol.add_atom(Atom::new("H", [1.0; 3]));
    let a3 = mol.add_atom(Atom::new("H", [2.0; 3]));
    mol.add_bond(a1, a2, Bond::default());
    mol.add_bond(a1, a3, Bond::default());
    // update bond type
    mol.add_bond(a3, a1, Bond::double());
    let bv = BondsView::new(&mol);
    {
        let b12 = &bv[(1, 2)];
        let b13 = &bv[(1, 3)];
    }

    // a list of tuple: (1, 2, Bond)
    let ijbs: Vec<_> = bv.collect();
    assert_eq!(2, ijbs.len());
    let (i, j, _) = ijbs[0];
    assert_eq!((1, 2), (i, j));
}

impl Molecule {
    pub fn view_bonds(&self) -> BondsView {
        BondsView::new(&self)
    }
}
// bonds view:1 ends here
