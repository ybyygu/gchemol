// [[file:~/Workspace/Programming/gchemol/gchemol.note::8842a219-a252-4367-bb8a-7a28b6bb8c2f][8842a219-a252-4367-bb8a-7a28b6bb8c2f]]
use Atom;
use Molecule;

use parser::{
    space_token,
    take_until_end_of_line,
    alphanumeric,
    digit_one_line,
    xyz_array,
};

named!(get_atom_from<&str, Atom>, do_parse!(
    // element symbol, "1" or "H"
    sym      : sp!(alphanumeric) >>
    position : sp!(xyz_array)    >>
    // ignore the remaining characters
    take_until_end_of_line       >>
    (
        Atom::new(sym, position)
    )
));

#[test]
fn test_format_xyz_atom() {
    let (_, x) = get_atom_from("C -11.4286 -1.3155  0.0000 ").unwrap();
    assert_eq!("C", x.symbol());
}

/// parse molecule from &str
named!(get_molecule_from<&str, Molecule>, do_parse!(
    natoms: sp!(digit_one_line)    >>
    title : take_until_end_of_line >>
    atoms : many0!(get_atom_from)  >>
    (
        {
            let mut mol = Molecule::new(title);
            if atoms.len() != natoms {
                eprintln!("the expected number of atoms is different: {}, {}", natoms, atoms.len());
            }

            for a in atoms {
                mol.add_atom(a);
            }
            mol
        }
    )
));

#[test]
fn test_formats_xyz_molecule() {
    let txt = "12
title
C -11.4286  1.7645  0.0000
C -10.0949  0.9945  0.0000
C -10.0949 -0.5455  0.0000
C -11.4286 -1.3155  0.0000
C -12.7623 -0.5455  0.0000
C -12.7623  0.9945  0.0000
H -11.4286  2.8545  0.0000
H -9.1509  1.5395  0.0000
H -9.1509 -1.0905  0.0000
H -11.4286 -2.4055  0.0000
H -13.7062 -1.0905  0.0000
H -13.7062  1.5395  0.0000";

    let (_, mol) = get_molecule_from(txt).unwrap();
    assert_eq!(12, mol.natoms());
}
// 8842a219-a252-4367-bb8a-7a28b6bb8c2f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c6258370-89a6-4eda-866c-41d60ef03e44][c6258370-89a6-4eda-866c-41d60ef03e44]]
use std::str;
use std::fs::File;
use std::io::prelude::*;

use nom::IResult;

use errors::*;
use formats::{
    ChemFileLike,
};

pub struct XYZFile();

impl ChemFileLike for XYZFile {
    fn ftype(&self) -> &str {
        "text/xyz"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".xyz"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule_from(chunk)
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        // meta information
        let mut lines = String::new();
        lines.push_str(&format!("{}\n", mol.natoms()));
        lines.push_str(&format!("{}\n", mol.title()));

        // coordinates
        for a in mol.atoms() {
            let p = a.position();
            let sym = a.symbol();
            let s = format!("{:6} {:-18.6}{:-18.6}{:-18.6}\n",
                            sym,
                            p[0], p[1], p[2]);
            lines.push_str(&s);
        }

        Ok(lines)
    }
}

#[test]
fn test_formats_xyz() {
    let file = XYZFile();
    // let x = file.parse("/home/ybyygu/Workspace/Projects/structure-prediction/nanoreactor/tests/deMon.xyz");
    let mols = file.parse("tests/files/xyz/multi.xyz").unwrap();
    assert_eq!(6, mols.len());
}
// c6258370-89a6-4eda-866c-41d60ef03e44 ends here
