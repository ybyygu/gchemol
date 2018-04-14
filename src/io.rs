// [[file:~/Workspace/Programming/gchemol/gchemol.note::891f59cf-3963-4dbe-a7d2-48279723b72e][891f59cf-3963-4dbe-a7d2-48279723b72e]]
//===============================================================================#
//   DESCRIPTION:  basic read & write support for molecular file
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-11 Wed 15:42>
//       UPDATED:  <2018-04-14 Sat 18:45>
//===============================================================================#
// 891f59cf-3963-4dbe-a7d2-48279723b72e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::0f52a1ef-c664-45a9-ab96-6d31741ae8c0][0f52a1ef-c664-45a9-ab96-6d31741ae8c0]]
use std::io::prelude::*;
use std::io::{BufWriter, BufReader};
use std::fs::File;

use errors::*;
use Point3D;
use Points;

/// Return content of text file in string
/// don't use if file is large
pub fn read_file(filename: &str) -> Result<String> {
    let mut buffer = String::new();

    let mut fp = File::open(filename).chain_err(|| format!("unable to open {} for reading", filename))?;
    fp.read_to_string(&mut buffer).chain_err(||format!("failed to read content: {}", filename))?;

    Ok(buffer)
}

/// write string content to an external file
pub fn write_file(content: String, filename: &str) -> Result<()> {
    let msg = format!("failed to create output file: {}", filename);
    let f = File::create(filename)
        .chain_err(|| msg)?;
    let mut writer = BufWriter::new(&f);

    let msg = format!("failed to write output file: {}", filename);
    writer.write_all(&content.as_bytes())
        .chain_err(|| msg)?;

    Ok(())
}

/// write a list of string without line ending characters to an external file
pub fn write_lines(lines: Vec<String>, filename: &str) -> Result<()> {
    let content = lines.join("\n");
    write_file(content, filename);

    Ok(())
}
// 0f52a1ef-c664-45a9-ab96-6d31741ae8c0 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::618eea82-a702-4d2f-a873-3807ead50d4b][618eea82-a702-4d2f-a873-3807ead50d4b]]
/// write coordinates in xyz format
pub fn write_as_xyz(symbols: &[&str], positions: &Points, filename: &str) -> Result<()>
{
    let mut lines = String::new();

    // meta information
    lines.push_str(format!("{}\n", positions.len()).as_str());
    lines.push_str("default title\n");

    // coordinates
    for i in 0..symbols.len() {
        let p = positions[i];
        let sym = symbols[i];
        let s = format!("{:6} {:-18.6}{:-18.6}{:-18.6}\n", sym, p[0], p[1], p[2]);
        lines.push_str(s.as_str());
    }

    // final empty line
    lines.push_str("\n");

    write_file(lines, filename);

    Ok(())
}
// 618eea82-a702-4d2f-a873-3807ead50d4b ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a1e0083b-88b2-46d1-ae09-3ad7165fdb9e][a1e0083b-88b2-46d1-ae09-3ad7165fdb9e]]
/// read essential molecular data from a xyz file
fn read_xyzfile(filename: &str) -> Result<(Vec<String>, Points)> {
    let text = read_file(filename)?;
    let mut lines: Vec<_> = text.lines().collect();

    let nlines = lines.len();
    let natoms: usize = lines[0].parse()
        .chain_err(|| format!("failed to parse int number from: {:?}", lines[0]))?;
    if natoms != nlines - 2 {
        eprintln!("the expected number of atoms is inconsistent with the xyz records.");
    }

    let mut positions = vec![];
    let mut symbols = vec![];
    for line in &lines[2..(natoms+2)] {
        let attrs: Vec<_> = line.split_whitespace().collect();
        let (symbol, position) = attrs.split_first().ok_or("encountering empty line")?;
        if position.len() != 3 {
            let msg = format!("informal xyz records: {}", line);
            bail!(msg);
        }

        symbols.push(symbol.to_string());

        let p: Vec<f64> = position.iter().map(|x| x.parse().unwrap()).collect();
        positions.push([p[0], p[1], p[2]]);

    }

    Ok((symbols, positions))
}

#[test]
#[ignore]
fn test_read_xyzfile() {
    let (symbols, positions) = read_xyzfile("/tmp/test.xyz").unwrap();
    println!("{:?}", symbols);
}
// a1e0083b-88b2-46d1-ae09-3ad7165fdb9e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::01bb2894-4548-42f3-9977-84164657219c][01bb2894-4548-42f3-9977-84164657219c]]
use std::collections::HashMap;

/// A quick and dirty way to get molecules from .mol2 file
///
/// Parameters
/// ---------
/// filename: the path to a tripos mol2 file
///
/// Return
/// ------
/// A list of molecules
fn from_mol2file(filename: &str) -> Result<Molecule> {
    let txt = read_file(filename)?;

    let mut lines = txt.lines();

    let mut jumping_to = |tag: &str| {
    };

    //
    // 1. get natoms and nbonds
    // -----------------------------------------------------------------------------
    // jumping
    let tag = "@<TRIPOS>MOLECULE";
    loop {
        if let Some(line) = lines.next() {
            if line.trim_left().contains(tag) {
                break;
            }
        } else {
            bail!("expected tag {} not found", tag);
        }
    }

    // skip one line and get the next line
    let mut natoms = 0;
    let mut nbonds = 0;
    if let (_, Some(line)) = (lines.next(), lines.next()) {
        let parts: Vec<_> = line.split_whitespace().collect();
        if parts.len() == 2 {
            natoms = parts[0].parse().chain_err(|| "cannot get natoms")?;
            nbonds = parts[1].parse().chain_err(|| "cannot get nbonds")?;
        } else {
            bail!("wrong line: {}", line);
        }
    } else {
        bail!("cannot read number of atoms and bonds in: {}", filename);
    }

    // 2. read in atoms
    let tag = "@<TRIPOS>ATOM";
    loop {
        if let Some(line) = lines.next() {
            if line.trim_left().contains(tag) {
                break;
            }
        } else {
            bail!("expected tag {} not found", tag);
        }
    }

    // 2. read in atoms
    let mut atoms = HashMap::new();
    let mut molecule = Molecule::new();
    // atom_id element_symbol, x, y, z
    for _ in 0..natoms {
        if let Some(line) = lines.next() {
            let parts: Vec<_> = line.split_whitespace().take(5).collect();
            let index = parts[0];
            let label = parts[1];
            // remove any numbers after element symbol: e.g. N39
            let symbol = label.trim_right_matches(char::is_numeric);
            let x: f64 = parts[2].parse().unwrap();
            let y: f64 = parts[3].parse().unwrap();
            let z: f64 = parts[4].parse().unwrap();
            let mut atom = Atom::new(symbol, [x, y, z]);
            atom.name = label.to_string();
            let a = molecule.add_atom(atom);
            atoms.insert(index, a);
        } else {
            bail!("incomplete atom records");
        }
    }

    // 3. read in bonds
    let tag = "@<TRIPOS>BOND";
    loop {
        if let Some(line) = lines.next() {
            if line.trim_left().contains(tag) {
                break;
            }
        } else {
            bail!("expected tag {} not found", tag);
        }
    }
    for _ in 0..nbonds {
        if let Some(line) = lines.next() {
            // ignore bond order attribute
            let parts: Vec<_> = line.split_whitespace()
                .skip(1)
                .take(2)
                .collect();
            let i = parts[0];
            let j = parts[1];
            let ai = atoms[i];
            let aj = atoms[j];
            molecule.add_bond(ai, aj);
        } else {
            bail!("incomplete bond records");
        }
    }

    Ok(molecule)
}

#[test]
fn test_from_mol2file() {
    // FIXME: the path
    let mol = from_mol2file("tests/files/alanine-gv.mol2").unwrap();
    assert_eq!(12, mol.natoms());
    assert_eq!(11, mol.nbonds());
}
// 01bb2894-4548-42f3-9977-84164657219c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::bd8ea4f1-b38f-4c81-a722-96cf78fb53a6][bd8ea4f1-b38f-4c81-a722-96cf78fb53a6]]
fn to_mol2file() {
    //
}
// bd8ea4f1-b38f-4c81-a722-96cf78fb53a6 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::448a8479-f0e8-412a-9e8d-83865581eb43][448a8479-f0e8-412a-9e8d-83865581eb43]]
use {
    Atom,
    Molecule,
};

impl Molecule {
    /// Construct molecule from external text file
    pub fn from_file<T: Into<String>>(filename: T) -> Self {
        let filename = filename.into();
        let (symbols, positions) = read_xyzfile(&filename).unwrap();
        let mut mol = Molecule::new();
        for i in 0..symbols.len() {
            let sym = &symbols[i];
            let pos = &positions[i];
            let atom = Atom::new(sym.clone(), pos.clone());
            mol.add_atom(atom);
        }

        mol
    }

    /// Save molecule to an external file
    pub fn to_file<T: Into<String>>(&self, filename: T) -> Result<()> {
        let filename = filename.into();

        let symbols: Vec<_> = self.symbols().collect();
        let positions: Vec<_> = self.positions().map(|v| *v).collect();

        write_as_xyz(&symbols, &positions, &filename);
        Ok(())
    }
}

#[test]
fn test_molecule_from_file() {
    let mol = Molecule::from_file("tests/data/c2h4.xyz");
    assert_eq!(6, mol.natoms());
}
// 448a8479-f0e8-412a-9e8d-83865581eb43 ends here
