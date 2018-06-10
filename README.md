gchemol is a graph-based chemical object library implemented in Rust programming
language.

[![Build Status](https://travis-ci.org/ybyygu/gchemol.svg?branch=master)](https://travis-ci.org/ybyygu/gchemol)
[![GPL3 licensed](https://img.shields.io/badge/license-GPL3-blue.svg)](./LICENSE)


# Features

-   Fast and safe.
-   Easy to deploy in server environment.
-   core graph data structure using [petgraph](https://github.com/bluss/petgraph)
-   read molecules in various formats using [nom](https://github.com/Geal/nom) parser combinators
-   linear algebra backed by [nalgebra](http://nalgebra.org/)
-   render molecule in user defined formats by templating with [handlebars](https://github.com/sunng87/handlebars-rust)


# Installing


## install rust

follow the official guide:

-   [Installation · The Rust Programming Language](https://www.rust-lang.org/en-US/install.html)

add gchemol dependency to your Cargo.toml:

    [dependencies]
    gchemol = {git = "https://github.com/ybyygu/gchemol"}


# Edit molecule


## atom

    use gchemol::Atom;
    
    // construct from element and position
    let a = Atom::new("C", [0.0, 0.0, 0.0]);
    let b = Atom::new("C", [1.2, 1.2, 1.2]);

convert from string:

    let a: Atom = "C 1.0 1.0 0.2"
        .parse()
        .expect("atom from string");

using the builder pattern

    let a = Atom::build()
        .symbol("Fe")
        .position([1.2, 1.0, 0.3])
        .finish();


## molecule

1.  Creating a molecule manually

        use gchemol::Molecule;
        
        let atoms = [
            Atom::new("C", [ 0.000000,   0.000000,  0.000000]),
            Atom::new("C", [ 0.000000,   0.000000,  1.089000]),
            Atom::new("C", [ 1.026719,   0.000000, -0.363000]),
            Atom::new("C", [-0.513360,  -0.889165, -0.363000]),
            Atom::new("C", [-0.513360,   0.889165, -0.363000])];
        
        // create a molecule named methane
        let mut mol = Molecule::new("methane");
        // add atoms in a loop
        for a in atoms {
            mol.add_atom(a);
        }

2.  Reading and writing molecules

    The format is automatically guessed from the filename.
    
        use gchemol::io;
        use gchemol::Molecule;
        
        // Read an xyz file and write to a Gaussian Input file.
        let mol = Molecule::from_file("path/to/file");
        mol.to_file("methane.gjf")
        
        // get the total number of atoms
        let na = mol.natoms();
        // get the total number of bonds
        let nb = mol.nbonds();
        
        // read multiple molecules (trajectory) from a chemical file
        // the number of atoms in different frame could be different
        let mols = io::read("path/to/trajectory.xyz");

3.  coordinates

        let mut positions = mol.positions();
        positions[0] = [1.2, 1.0, 0.1];
        mol.set_positions(positions);
    
    to be continued &#x2026;


## bonds

1.  bonding connectivity

    TBD


## lattice

build a periodic structure

    use gchemol::Lattice;
    
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
    
    // set periodic boundary conditions for molecule
    mol.set_lattice(lat);
    
    // remove periodic boundary conditions
    mol.unbuild_crystal();


# Read/write chemical files

TBD


# Templating

For syntax details: [sunng87/handlebars-rust: Rust templating with Handlebars](https://github.com/sunng87/handlebars-rust)

molecule title

    {{molecule_title}}

total number of atoms in molecule

    {{molecule.number_of_atoms}}

total number of bonds in molecule

    {{molecule.number_of_bonds}}

loop over atoms

    {{#each molecule.atoms as |a| ~}}
    {{a.symbol}} {{a.x}} {{a.y}} {{a.z}}
    {{/each~}}

atom properties

    {{atom.symbol}}

cartesian coordinates

    {{atom.x}} {{atom.y}} {{atom.z}}

fractional coordinates

    {{atom.fx}} {{atom.fy}} {{atom.fz}}

unit cell parameters

    {{molecule.unit_cell.a}}
    {{molecule.unit_cell.b}}
    {{molecule.unit_cell.c}}
    {{molecule.unit_cell.alpha}}
    {{molecule.unit_cell.beta}}
    {{molecule.unit_cell.gamma}}
    {{molecule.unit_cell.va}}
    {{molecule.unit_cell.vb}}
    {{molecule.unit_cell.vc}}

format string or number:

    {{format 12.16 width=12 prec=6 align="right"}}

alignment control:

    {{align="left"}}
    {{align="right"}}
    {{align="center"}}


# Related projects

-   lumol
-   ase
-   pymatgen


# References

-   [Lessons from sixteen years of molecular simulation in Python | Konrad Hinsen's Blog](https://khinsen.wordpress.com/2013/04/10/lessons-from-sixteen-years-of-molecular-simulation-in-python/)

