// build lattice

use approx::*;

#[test]
fn test_crystal() {
    use gchemol::{
        Molecule,
        Lattice,
    };

    // build lattice from lattice parameters
    let _ = Lattice::from_params
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
    use gchemol:: {
        AtomIndex,
        Molecule,
        prelude::*,
    };

    use gchemol::io;

    let mut mols = io::read("tests/files/cif/MS-MOR.cif").expect("cif test file");
    let mol = &mut mols[0];
    let d = mol.distance(AtomIndex::new(0), AtomIndex::new(12)).expect("distance between 0 and 12");
    assert_relative_eq!(12.6753, d, epsilon=1e-4);

    // remove periodic bound
    mol.unbuild_crystal();
    let d = mol.distance(AtomIndex::new(0), AtomIndex::new(12)).expect("distance between 0 and 12");
    assert_relative_eq!(16.203993, d, epsilon=1e-4);

    // distance matrix
    let mol = Molecule::from_file("tests/files/cif/quinone.cif").expect("mol2 test file");

    let dm = mol.distance_matrix();
    let expected = [0.        , 1.46623609, 1.46701952, 1.21833857, 4.73857325,
                    4.76062156, 4.00902774, 4.20621928, 2.82936072, 2.41889587,
                    2.42331482, 4.04768958, 3.69004223, 3.4732701 , 4.16339654,
                    3.34305314];
    for i in 0..expected.len() {
        assert_relative_eq!(expected[i], dm[i], epsilon=1e-4);
    }
}

// supercell

#[test]
fn test_supercell() -> Result<()> {
    use gchemol::io::prelude::*;
    use gchemol::Supercell;

    let path = "/tmp/a.cif";
    let images = gchemol::io::read(path)?;

    let mut images_modified = vec![];
    for image in images {
        let mol = Supercell::new()
            .with_range_a(-1, 1)
            .with_range_b(-1, 1)
            .build(&image);
        images_modified.push(mol);
    }
    gchemol::io::write("/tmp/all.cif", &images_modified)?;

    Ok(())
}
