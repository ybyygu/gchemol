// demonstrate how to build molecule manually

// [[file:~/Workspace/Programming/gchemol/gchemol.note::*demonstrate%20how%20to%20build%20molecule%20manually][demonstrate how to build molecule manually:1]]
extern crate gchemol;
#[macro_use] extern crate approx;
// demonstrate how to build molecule manually:1 ends here

// construct atom with element and position

// [[file:~/Workspace/Programming/gchemol/gchemol.note::*construct%20atom%20with%20element%20and%20position][construct atom with element and position:1]]
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
        .momentum(0.0, 5.0, 6.0)
        .finish();

    assert_eq!(a.symbol(), "Fe");
    assert_eq!(a.position(), [1.2, 1.0, 0.3]);
    assert_eq!([0.0, 5.0, 6.0], a.momentum());
}
// construct atom with element and position:1 ends here

// build molecule from atoms

// [[file:~/Workspace/Programming/gchemol/gchemol.note::*build%20molecule%20from%20atoms][build molecule from atoms:1]]
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
    assert_eq!(5, mol.natoms());

    // update positions
    let positions = [[-0.90203687,  0.62555259,  0.0081889 ],
                     [-0.54538244, -0.38325741,  0.0081889 ],
                     [-0.54536403,  1.12995078,  0.8654626 ],
                     [-0.54536403,  1.12995078, -0.8654626 ],
                     [-1.97203687,  0.62556577,  0.0081889 ]];
    mol.set_positions(&positions).expect("atom positions");

    // udpate symbols
    // accept array
    let symbols = ["C", "H", "H", "H", "H"];
    mol.set_symbols(&symbols);
    // accept vec
    let symbols = vec!["C", "H", "H", "H", "H"];
    mol.set_symbols(&symbols);
    for (s1, s2) in mol.symbols().into_iter().zip(symbols) {
        assert_eq!(s1, s2);
    }
}
// build molecule from atoms:1 ends here

// properties
// Set arbitrary properties for atom.

// [[file:~/Workspace/Programming/gchemol/gchemol.note::*properties][properties:1]]
#[test]
fn test_molecule_properties() {
    use gchemol::Atom;

    // create Atom object
    let mut a = Atom::default();

    // store a key with an integer value
    let value = 12;
    a.properties.store("secret", value);
    let unpacked: isize = a.properties.load("secret").expect("secret property");
    assert_eq!(value, unpacked);

    // store a key with a list of integers
    let value = [1, 2, 3];
    a.properties.store("secret", value);

    // test if the property if exists
    assert!(a.properties.contains_key("secret"));
    assert!(!a.properties.contains_key("blank"));

    // unpack the property
    let unpacked: Vec<isize> = a.properties.load("secret").expect("secret property");
    for i in 0..value.len() {
        assert_eq!(value[i], unpacked[i]);
    }
}
// properties:1 ends here

// query bonded atoms

// [[file:~/Workspace/Programming/gchemol/gchemol.note::*query%20bonded%20atoms][query bonded atoms:1]]
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
// query bonded atoms:1 ends here

// molecule center

// [[file:~/Workspace/Programming/gchemol/gchemol.note::*molecule%20center][molecule center:1]]
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
// molecule center:1 ends here
