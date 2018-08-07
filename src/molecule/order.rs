// [[file:~/Workspace/Programming/gchemol/gchemol.note::a762197d-df95-433c-8499-8148d0241a9f][a762197d-df95-433c-8499-8148d0241a9f]]
use super::*;
use std::cmp::Reverse;

/// reversely sort atoms by atom numbers, but preserves their insertion order
/// for atoms with the same element type
fn rev_sort_atoms_by_element(mol: &Molecule) -> Vec<(usize, AtomIndex)> {
    let numbers = mol.numbers();
    let sites = mol.sites();

    let mut pairs: Vec<_> = numbers.into_iter().zip(sites.into_iter()).collect();
    pairs.sort_by_key(|&ns| (ns.0, Reverse(ns.1)));
    // hydrogen last
    pairs.reverse();

    pairs
}

impl Molecule {
    // FIXME: keep or remove?
    /// Update atom indices
    pub fn reorder(&mut self) {
        let ns: Vec<_> = self.graph.node_indices().collect();
        for (i, &n) in ns.iter().enumerate() {
            let atom = &mut self.graph[n];
            atom.index = n;
        }
    }

    /// Return a new molecule with all atoms sorted by element number (hydrogen last)
    pub fn sorted(&self) -> Molecule {
        let title = format!("sorted {}", self.title());
        let mut mol = Molecule::new(&title);

        let pairs = rev_sort_atoms_by_element(&self);

        let mut mapping = HashMap::with_capacity(pairs.len());
        for (_, n) in pairs {
            let a = self.get_atom(n).expect("sorted: pairs");
            let new_n = mol.add_atom(a.clone());
            mapping.insert(n, new_n);
        }

        for b in self.bonds() {
            let e = b.index();
            let (o1, o2) = self.partners(&e).expect("sorted: bonds");
            let (n1, n2) = (mapping[&o1], mapping[&o2]);
            mol.add_bond(n1, n2, b.clone());
        }

        // lattice
        if let Some(lat) = self.lattice {
            mol.set_lattice(lat.clone());
        }

        // properties
        mol.properties = self.properties.clone();

        mol
    }
}

#[test]
fn test_molecule_sorted() {
    let mol = Molecule::from_file("tests/files/mol2/alanine-gv.mol2").expect("mol2 for reorder");
    let mol_sorted = mol.sorted();
    assert_eq!(mol.natoms(), mol_sorted.natoms());
    assert_eq!(mol.nbonds(), mol_sorted.nbonds());
    let numbers = mol_sorted.numbers();
    assert_eq!(8, numbers[0]);
    assert_eq!(7, numbers[1]);
    assert_eq!(6, numbers[2]);
    assert_eq!(6, numbers[3]);
}
// a762197d-df95-433c-8499-8148d0241a9f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b5512aff-1510-42cf-9b1d-7487485a282d][b5512aff-1510-42cf-9b1d-7487485a282d]]
// take from: https://rosettacode.org/wiki/Permutations#Rust
pub fn permutations(size: usize) -> Permutations {
    Permutations { idxs: (0..size).collect(), swaps: vec![0; size], i: 0 }
}

pub struct Permutations {
    idxs: Vec<usize>,
    swaps: Vec<usize>,
    i: usize,
}

impl Iterator for Permutations {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.i > 0 {
            loop {
                if self.i >= self.swaps.len() { return None; }
                if self.swaps[self.i] < self.i { break; }
                self.swaps[self.i] = 0;
                self.i += 1;
            }
            self.idxs.swap(self.i, (self.i & 1) * self.swaps[self.i]);
            self.swaps[self.i] += 1;
        }
        self.i = 1;
        Some(self.idxs.clone())
    }
}

#[test]
fn test_permutation() {
    let perms: Vec<_> = permutations(3).collect();
    assert_eq!(perms, vec![
        vec![0, 1, 2],
        vec![1, 0, 2],
        vec![2, 0, 1],
        vec![0, 2, 1],
        vec![1, 2, 0],
        vec![2, 1, 0],
    ]);
}
// b5512aff-1510-42cf-9b1d-7487485a282d ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::e06ad932-0eef-4d99-beac-c43c4e83bc63][e06ad932-0eef-4d99-beac-c43c4e83bc63]]
impl Molecule {
    /// Test if other molecule suitable for matching.
    pub fn matchable(&self, other: &Molecule) -> bool {
        // The two molecules must contain the same numbers of atoms of each
        // element type.
        if self.natoms() == other.natoms() {
            let reduced_symbols1 = self.reduced_symbols();
            let reduced_symbols2 = other.reduced_symbols();
            reduced_symbols1 == reduced_symbols2
        } else {
            false
        }
    }
}

#[test]
fn test_molecule_match() {
    let mol_ref = Molecule::from_file("tests/files/alignment/reference.mol2").expect("mol2 reference");
    let mol_can = Molecule::from_file("tests/files/alignment/candidate-ro.mol2").expect("mol2 candidate");
    assert!(mol_ref.matchable(&mol_can));
}
// e06ad932-0eef-4d99-beac-c43c4e83bc63 ends here
