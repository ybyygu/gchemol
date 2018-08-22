// build lattice

extern crate gchemol;
#[macro_use] extern crate approx;

use gchemol::Lattice;

fn test_crystal() {
    use gchemol::{
        Molecule,
        Lattice,
    };

    // build lattice from lattice parameters
    let lattice = Lattice::from_params
        (
            3.84,                   // a
            3.84,                   // b
            3.84,                   // c
            120.,                   // alpha
            90.,                    // beta
            60.,                    // gamma
        );

    // build lattice from cell vectors
    let lat = Lattice::new([
        [ 15.3643,   0.    ,   0.    ], // vector a
        [  4.5807,  15.5026,   0.    ], // vector b
        [  0.    ,   0.    ,  17.4858], // vector c
    ]);

    // create an empty molecule
    let mut mol = Molecule::new("empty");

    // set periodic boundary conditions for molecule
    mol.set_lattice(lat);

    // remove periodic boundary conditions
    mol.unbuild_crystal();
}

// distance using mic

#[test]
fn test_molecule_pbc_distance() {
    use gchemol::AtomIndex;
    use gchemol::io;

    let mut mols = io::read("tests/files/cif/MS-MOR.cif")
        .expect("structure from cif file");
    let mut mol = &mut mols[0];
    let d = mol.distance(AtomIndex::new(0), AtomIndex::new(12))
        .expect("distance between 0 and 12");
    assert_relative_eq!(12.6753, d, epsilon=1e-4);

    // remove periodic bound
    mol.unbuild_crystal();
    let d = mol.distance(AtomIndex::new(0), AtomIndex::new(12))
        .expect("distance between 0 and 12");
    assert_relative_eq!(16.203993, d, epsilon=1e-4);
}
