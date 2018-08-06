// [[file:~/Workspace/Programming/gchemol/gchemol.note::a762197d-df95-433c-8499-8148d0241a9f][a762197d-df95-433c-8499-8148d0241a9f]]
use super::*;

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

        let numbers = self.numbers();
        let sites = self.sites();
        let mut pairs: Vec<_> = numbers.iter().zip(sites.iter()).collect();
        pairs.sort();
        pairs.reverse();

        let mut mapping = HashMap::with_capacity(pairs.len());
        for (_, &n) in pairs {
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
