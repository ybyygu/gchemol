// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*atom/molecule][atom/molecule:1]]
use super::*;

named!(get_atom_from<&str, Atom>, do_parse!(
    // element symbol, "1" or "H"
    sym      : sp!(alt!(alpha|digit)) >>
    position : sp!(xyz_array)         >>
    // ignore the remaining characters
               read_until_eol         >>
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

named!(get_molecule_pxyz<&str, Molecule>, do_parse!(
    atoms: many1!(get_atom_from)              >>
           // match the molecule record separator
           alt!(blank_line | tag!(MAGIC_EOF)) >>
    (
        {
            let mut mol = Molecule::new("untitled molecule");
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
H -13.7062  1.5395  0.0000
";
    let (_, mol) = get_molecule(txt).unwrap();
    assert_eq!(12, mol.natoms());
}
// atom/molecule:1 ends here

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*chemfile][chemfile:1]]
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


/// plain xyz coordinates with atom symbols
#[derive(Debug, Clone)]
pub struct PlainXYZFile();

impl ChemFileLike for PlainXYZFile {
    /// possible file extensions
    fn extensions(&self) -> Vec<&str> {
        [".coord", ".pxyz", ".coords"].to_vec()
    }

    fn ftype(&self) -> &str {
        "text/pxyz"
    }

    /// Parse a single molecule from file `filename`
    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule_pxyz(chunk)
    }

    /// Return a string representation of molecule
    /// Multiple molecules will be separated by a blank line
    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        let mut lines = String::new();

        for a in mol.atoms() {
            lines.push_str(format!("{}\n", a.to_string()).as_ref());
        }

        // append a blank line as a separator between multiple molecules
        lines.push_str("\n");

        Ok(lines)
    }
}

#[test]
fn test_formats_plain_xyz() {
    let filename = "tests/files/xyz/c2h4.pxyz";
    let file = PlainXYZFile();
    assert!(file.parsable(filename));
    let mols = file.parse(filename).unwrap();
    assert_eq!(1, mols.len());
    assert_eq!(6, mols[0].natoms());

    // parse multiple molecules
    let mols = file.parse("tests/files/xyz/multi.pxyz").expect("multi xyz");
    assert_eq!(6, mols.len());

    let natoms_expected = vec![16, 10, 16, 16, 16, 13];
    let natoms: Vec<_> = mols.iter().map(|m| m.natoms()).collect();
    assert_eq!(natoms_expected, natoms);
}

#[test]
fn test_formats_xyz() {
    let file = XYZFile();
    let mols = file.parse("tests/files/xyz/c2h4.xyz").expect("c2h4 xyz");
    assert_eq!(1, mols.len());
    assert_eq!(6, mols[0].natoms());

    // parse multiple molecules
    let mols = file.parse("tests/files/xyz/multi.xyz").expect("multi xyz");
    assert_eq!(6, mols.len());

    let natoms_expected = vec![16, 10, 16, 16, 16, 13];
    let natoms: Vec<_> = mols.iter().map(|m| m.natoms()).collect();
    assert_eq!(natoms_expected, natoms);
}
// chemfile:1 ends here
