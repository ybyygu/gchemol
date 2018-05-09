// [[file:~/Workspace/Programming/gchemol/gchemol.note::ea19d54a-dbe6-4367-90ea-5a0465018219][ea19d54a-dbe6-4367-90ea-5a0465018219]]
named!(cell_params<&str, (f64, f64, f64, f64, f64, f64)>, permutation!(
    preceded!(ws!(tag!("_cell_length_a")), ws!(double_s)),
    preceded!(ws!(tag!("_cell_length_b")), ws!(double_s)),
    preceded!(ws!(tag!("_cell_length_c")), ws!(double_s)),
    preceded!(ws!(tag!("_cell_angle_alpha")), ws!(double_s)),
    preceded!(ws!(tag!("_cell_angle_beta")), ws!(double_s)),
    preceded!(ws!(tag!("_cell_angle_gamma")), ws!(double_s))
));

#[test]
fn test_cif_cell_loop() {
    let lines = "_cell_length_a                    18.0940
_cell_length_c                    7.5240
_cell_length_b                    20.5160
_cell_angle_alpha                 90.0000
_cell_angle_beta                  90.0000
_cell_angle_gamma                 90.0000
";

    let (_, x) = cell_params(lines).unwrap();
    assert_relative_eq!(20.5160, x.1, epsilon=1e-3);
}
// ea19d54a-dbe6-4367-90ea-5a0465018219 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::6dcc25a3-6738-497a-9317-df051c7afa74][6dcc25a3-6738-497a-9317-df051c7afa74]]
use Atom;

use std::collections::HashMap;

use parser::{
    alpha,
    alphanumeric,
    double_s,
    take_until_end_of_line,
    space_token,
    end_of_line,
    not_space,
};

use nom::IResult;

named!(atom_site_header<&str, &str>, preceded!(
    tag!("_atom_site_"),
    not_space
));

#[test]
fn test_cif_atom_site_header() {
    let (_, x) = atom_site_header("_atom_site_type_symbol").unwrap();
    assert_eq!("type_symbol", x);
}

// ["label", "type_symbol", "fract_x", "fract_y", "fract_z", "U_iso_or_equiv", "adp_type", "occupancy"]
named!(atom_site_headers<&str, Vec<&str>>, do_parse!(
    ws!(tag!("loop_")) >>
    headers: many1!(ws!(atom_site_header)) >>
    (
        headers
    )
));

// How to deal with IResult function
// since I need handle the variable headers
fn atom_site_loops<'a>(input: &'a str) -> IResult<&'a str, Vec<Atom>> {
    // FIXME: handle possible errors
    let (rest, headers) = atom_site_headers(input).unwrap();

    assert!(headers.len() >= 4, "not enough columns in atom site loop");

    // column header loopup table
    // Example
    // -------
    //   0        1         2       3      4            5          6         7
    // label type_symbol fract_x fract_y fract_z U_iso_or_equiv adp_type occupancy
    let table: HashMap<_, _> = HashMap::from(headers.iter().zip(0..).collect());
    let ifx = *table.get(&"fract_x").expect("fract x col");
    let ify = *table.get(&"fract_y").expect("fract y col");
    let ifz = *table.get(&"fract_z").expect("fract z col");
    // column index to atom label
    let ilbl = *table.get(&"label").expect("atom label col");
    // TODO: column index to element symbol, which is optional
    let isym = *table.get(&"type_symbol").expect("atom symbol col");

    do_parse!(
        rest,
        rows: many1!(ws!(count!(sp!(not_space), headers.len()))) >>
        (
            {
                let mut atoms = vec![];
                for row in rows {
                    let  fx: f64 = row[ifx].parse().expect("cif fract x");
                    let  fy: f64 = row[ify].parse().expect("cif fract y");
                    let  fz: f64 = row[ifz].parse().expect("cif fract z");
                    let lbl = row[ilbl];
                    let sym = row[isym];
                    // TODO: assign atom label
                    let a = Atom::build()
                        .symbol(sym)
                        .position(fx, fy, fz)
                        .finish();
                    atoms.push(a);
                }

                atoms
            }
        )
    )
}

#[test]
fn test_cif_atoms_loop() {
    let lines = "loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy
Si1    Si    0.30070   0.07240   0.04120   0.00000  Uiso   1.00
Si2    Si    0.30370   0.30880   0.04610   0.00000  Uiso   1.00
O3     O     0.12430   0.41700   0.42870   0.00000  Uiso   1.00
O4     O     0.12260   0.19540   0.42540   0.00000  Uiso   1.00
O5     O     0.23620   0.12240   0.98650   0.00000  Uiso   1.00
Si6    Si    0.80070   0.57240   0.04120   0.00000  Uiso   1.00
Si7    Si    0.80370   0.80880   0.04610   0.00000  Uiso   1.00
O8     O     0.62430   0.91700   0.42870   0.00000  Uiso   1.00
O9     O     0.62260   0.69540   0.42540   0.00000  Uiso   1.00
O10    O     0.73620   0.62240   0.98650   0.00000  Uiso   1.00
Si11   Si    0.69930   0.92760   0.54120   0.00000  Uiso   1.00
Si12   Si    0.69630   0.69120   0.54610   0.00000  Uiso   1.00
";
    let x = atom_site_loops(lines);
    println!("{:?}", x);
}
// 6dcc25a3-6738-497a-9317-df051c7afa74 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::383cfb7d-0863-4efa-b969-4a6cbf7f3ad9][383cfb7d-0863-4efa-b969-4a6cbf7f3ad9]]
named!(geom_bond_header<&str, &str>, preceded!(
    alt!(
        tag!("_geom_bond_") |
        tag!("_ccdc_geom_bond_")
    ),
    not_space
));

named!(geom_bond_headers<&str, Vec<&str>>, preceded!(
    ws!(tag!("loop_")),
    many1!(ws!(geom_bond_header))
));

#[test]
fn test_cif_geom_bond_header() {
    let (_, x) = geom_bond_headers("loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type").unwrap();
    assert_eq!(5, x.len());
}
// 383cfb7d-0863-4efa-b969-4a6cbf7f3ad9 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c7482893-b288-449a-82ba-387c85f3e55c][c7482893-b288-449a-82ba-387c85f3e55c]]
use lattice::Lattice;

named!(get_molecule_from<&str, Molecule>, do_parse!(
    name   : take_until_end_of_line >>
             take_until1!("_cell_") >>
    params : cell_params            >>
             take_until1!("loop_\n_atom_") >>
    atoms  : atom_site_loops >>
    (
        {
            let mut mol = Molecule::new(name.trim());
            let lat = Lattice::from_params(
                params.0,
                params.1,
                params.2,
                params.3,
                params.4,
                params.5,
            );

            for mut a in atoms {
                let p = lat.to_cart(a.position());
                a.set_position(p);
                mol.add_atom(a);
            }

            mol.set_lattice(lat);

            mol
        }
    )
));

#[test]
fn test_cif_molecule() {
    let lines = " data_LTL
#**************************************************************************
#
# CIF taken from the IZA-SC Database of Zeolite Structures
# Ch. Baerlocher and L.B. McCusker
# Database of Zeolite Structures: http://www.iza-structure.org/databases/
#
# The atom coordinates and the cell parameters were optimized with DLS76
# assuming a pure SiO2 composition.
#
#**************************************************************************

_cell_length_a                  18.12600
_cell_length_b                  18.12600
_cell_length_c                   7.56700
_cell_angle_alpha               90.00000
_cell_angle_beta                90.00000
_cell_angle_gamma              120.00000

_symmetry_space_group_name_H-M     'P 6/m m m'
_symmetry_Int_Tables_number         191
_symmetry_cell_setting             hexagonal

loop_
_symmetry_equiv_pos_as_xyz
'+x,+y,+z'
'-y,+x-y,+z'
'-x+y,-x,+z'
'-x,-y,+z'
'+y,-x+y,+z'
'+x-y,+x,+z'
'-y,-x,+z'
'-x+y,+y,+z'
'+x,+x-y,+z'
'+y,+x,+z'
'+x-y,-y,+z'
'-x,-x+y,+z'
'-x,-y,-z'
'+y,-x+y,-z'
'+x-y,+x,-z'
'+x,+y,-z'
'-y,+x-y,-z'
'-x+y,-x,-z'
'+y,+x,-z'
'+x-y,-y,-z'
'-x,-x+y,-z'
'-y,-x,-z'
'-x+y,+y,-z'
'+x,+x-y,-z'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
    O1    O     0.2645    0.5289    0.2231
    O2    O     0.1099    0.4162    0.3263
    O3    O     0.1484    0.5742    0.2620
    O4    O     0.1365    0.4736    0.0000
    O5    O     0.0000    0.2797    0.5000
    O6    O     0.1628    0.3256    0.5000
    T1    Si    0.1648    0.4982    0.2030
    T2    Si    0.0959    0.3594    0.5000
";

    let (_, x) = get_molecule_from(lines).unwrap();
    assert_eq!(8, x.natoms());
}
// c7482893-b288-449a-82ba-387c85f3e55c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::246e575d-8f3b-419d-b007-febcd6e45991][246e575d-8f3b-419d-b007-febcd6e45991]]
use Molecule;
use super::ChemFileLike;
use errors::*;
use io;

pub struct CifFile();

impl ChemFileLike for CifFile {
    fn ftype(&self) -> &str {
        "text/cif"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".cif"]
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        unimplemented!()
    }

    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        let txt = io::read_file(filename).chain_err(|| "failed to open file")?;
        let (_, mol) = get_molecule_from(&txt).unwrap();

        Ok(vec![mol])
    }
}

#[test]
#[ignore]
fn test_formats_cif() {
    let mols = io::read("tests/files/cif/MS-MOR.cif");
}
// 246e575d-8f3b-419d-b007-febcd6e45991 ends here
