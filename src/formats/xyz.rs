// [[file:~/Workspace/Programming/gchemol/gchemol.note::8842a219-a252-4367-bb8a-7a28b6bb8c2f][8842a219-a252-4367-bb8a-7a28b6bb8c2f]]
use super::*;

named!(get_atom_from<&str, Atom>, do_parse!(
    // element symbol, "1" or "H"
    sym      : sp!(alt!(alpha|digit))      >>
    position : sp!(xyz_array)         >>
    // ignore the remaining characters
               read_until_eol >>
    (
        Atom::new(sym, position)
    )
));

#[test]
fn test_formats_xyz_atom() {
    let (_, x) = get_atom_from("C -11.4286 -1.3155  0.0000\n").unwrap();
    assert_eq!("C", x.symbol());
    let (_, x) = get_atom_from("6 -11.4286 -1.3155  0.0000 \n").unwrap();
    assert_eq!("C", x.symbol());
}

fn get_molecule(input: &str) -> IResult<&str, Molecule> {
    // 1. read number of atoms
    let (rest, n) = terminated!(input, sp!(unsigned_digit), line_ending)?;
    // 2. get molecule title
    let (rest, t) = read_until_eol(rest)?;
    let title = t.trim();
    // 3. collect atom records
    let (rest, atoms) = many_m_n!(rest, n, n, get_atom_from)?;
    // 4. construct molecule
    let mut mol = Molecule::new(title);
    for a in atoms {
        mol.add_atom(a);
    }

    Ok((rest, mol))
}

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
H -13.7062  1.5395  0.0000
";
    let (_, mol) = get_molecule(txt).unwrap();
    assert_eq!(12, mol.natoms());
}
// 8842a219-a252-4367-bb8a-7a28b6bb8c2f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c6258370-89a6-4eda-866c-41d60ef03e44][c6258370-89a6-4eda-866c-41d60ef03e44]]
use std::str;
use std::fs::File;
use std::io::prelude::*;

pub struct XYZFile();

impl ChemFileLike for XYZFile {
    fn ftype(&self) -> &str {
        "text/xyz"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".xyz"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule(chunk)
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
    let mols = file.parse("tests/files/xyz/multi.xyz").expect("multi xyz");
    assert_eq!(6, mols.len());
}


/// plain xyz coordinates with atom symbols
#[derive(Debug)]
pub struct PlainXYZFile();

impl ChemFileLike for PlainXYZFile {
    /// possible file extensions
    fn extensions(&self) -> Vec<&str> {
        [".coord"].to_vec()
    }

    fn ftype(&self) -> &str {
        "text/coord"
    }

    /// parse molecules from file `filename`
    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        let txt = io::read_file(filename)?;
        let mut mol = Molecule::new("from plain coordinates");
        for line in txt.lines() {
            let line = line.trim();
            if line.is_empty() {
                break;
            }
            let a: Atom = line.parse()?;
            mol.add_atom(a);
        }

        Ok(vec![mol])
    }

    /// Return a string representation of the last molecule in the list
    /// Return empty string if no molecule found
    fn format(&self, mols: &Vec<Molecule>) -> Result<String> {
        let mut lines = String::new();

        if let Some(mol) = mols.last() {
            for a in mol.atoms() {
                lines.push_str(format!("{}\n", a.to_string()).as_ref());
            }
        }

        Ok(lines)
    }
}

#[test]
fn test_formats_plainxyz() {
    let filename = "tests/files/plain-coords/test.coord";
    let file = PlainXYZFile();
    assert!(file.parsable(filename));
    let mols = file.parse(filename).unwrap();
    assert_eq!(1, mols.len());
    assert_eq!(12, mols[0].natoms());
}
// c6258370-89a6-4eda-866c-41d60ef03e44 ends here
