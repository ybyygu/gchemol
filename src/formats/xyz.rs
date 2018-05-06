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
use std::io::{self, BufReader};
use std::io::prelude::*;

use errors::*;
use formats::{
    ChemFileLike,
};

use nom::IResult;

const BUF_SIZE: usize = 8 * 1024;
fn read_large_file(filename: &str) -> Result<Vec<Molecule>> {
    let fp = File::open(filename).chain_err(|| "failed")?;
    let mut reader = BufReader::with_capacity(BUF_SIZE, fp);

    let mut mols: Vec<Molecule> = vec![];
    let mut remained = String::new();
    let mut chunk = String::new();
    'out: loop {
        let length = {
            let buffer = reader.fill_buf().chain_err(|| "file buffer reading error")?;
            let new = str::from_utf8(&buffer).unwrap().to_string();
            // chunk = buffer + remained
            chunk.clear();
            chunk.push_str(&remained);
            chunk.push_str(&new);
            loop {
                // println!("inner {:?}", j);
                // fill chunk with remained data
                match get_molecule_from(&chunk) {
                    IResult::Error(err) => {
                        eprintln!("{:?}", err);
                        eprintln!("{:}", chunk);
                        break 'out;
                    },
                    IResult::Done(r, mol) => {
                        // println!("got mol with {:?} atoms", mol.natoms());
                        mols.push(mol);
                        remained = String::from(r);
                    },
                    IResult::Incomplete(i) => {
                        // eprintln!("need data: {:?}", i);
                        remained = chunk.clone();
                        break
                    },
                }
                // clear chunk
                chunk.clear();
                chunk.push_str(&remained);
            }
            buffer.len()
        };

        if length == 0 { break; }
        reader.consume(length);
    }

    Ok(mols)
}

struct XYZFile();

impl ChemFileLike for XYZFile {
    fn ftype(&self) -> &str {
        "text/xyz"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".xyz"]
    }

    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        read_large_file(filename)
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
