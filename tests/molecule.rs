// demonstrate how to build molecule manually

extern crate gchemol;
#[macro_use] extern crate approx;

// construct atom with element and position

#[test]
fn test_new_atom() {
    use gchemol::Atom;

    let a = Atom::new("Fe", [0.0, 0.0, 0.0]);
    let b = Atom::new("C", [1.2, 1.2, 1.2]);
    assert_eq!(a.symbol(), "Fe");
    assert_eq!(b.position(), [1.2, 1.2, 1.2]);

    // or simply convert from a string:
    let a: Atom = "Fe 1.0 1.0 0.2"
        .parse()
        .expect("atom from string");
    assert_eq!(a.symbol(), "Fe");
    assert_eq!(a.position(), [1.0, 1.0, 0.2]);

    // set more attributes using the builder pattern
    let a = Atom::build()
        .symbol("Fe")
        .position(1.2, 1.0, 0.3)
        .finish();

    assert_eq!(a.symbol(), "Fe");
    assert_eq!(a.position(), [1.2, 1.0, 0.3]);
}

// build molecule from atoms

#[test]
fn test_new_molecule() {
    use gchemol::{
        Atom,
        Molecule,
    };

    let atoms = vec![
        Atom::new("C", [ 0.000000,   0.000000,  0.000000]),
        Atom::new("C", [ 0.000000,   0.000000,  1.089000]),
        Atom::new("C", [ 1.026719,   0.000000, -0.363000]),
        Atom::new("C", [-0.513360,  -0.889165, -0.363000]),
        Atom::new("C", [-0.513360,   0.889165, -0.363000])];

    // create a molecule named as methane
    let mut mol = Molecule::new("methane");
    assert_eq!(mol.name, "methane");

    // add build molecule using for loop
    for a in atoms {
        mol.add_atom(a);
    }
    assert_eq!(5, mol.natoms())
}

// query bonded atoms

#[test]
fn test_molecule_neighbors() {
    use gchemol::Molecule;
    use gchemol::prelude::FromFile;

    let mol = Molecule::from_file("tests/files/mol2/alanine-gv.mol2").unwrap();
    assert_eq!(12, mol.natoms());

    let atoms = mol.view_atoms();
    let a3 = &atoms[3];

    let ns = a3.neighbors(&mol);
    assert_eq!(4, ns.len());

    let bonds = mol.view_bonds();
    let b35 = &bonds[(3, 5)];

    let (p1, p2) = b35.partners(&mol).unwrap();
}

// molecule center

#[test]
fn test_molecule_center() {
    use gchemol::Molecule;
    use gchemol::prelude::FromFile;

    let mol = Molecule::from_file("tests/files/mol2/alanine-gv.mol2").expect("mol2 gv");
    let pc = mol.center_of_geometry();
    let pe = [-2.31413333, -1.24455833,  0.41005833];

    for i in 0..3 {
        assert_relative_eq!(pe[i], pc[i], epsilon=1e-4);
    }
}
