// header
// #+name: fb7c688e-38d6-455b-a8f0-ae54c563f3cf

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::fb7c688e-38d6-455b-a8f0-ae54c563f3cf][fb7c688e-38d6-455b-a8f0-ae54c563f3cf]]
// gchemol parses the following record types in a PDB file:
//
// CRYST
// ATOM & HETATM
// TER
// END
// CONECT
// fb7c688e-38d6-455b-a8f0-ae54c563f3cf ends here

// base

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*base][base:1]]
use super::*;
// base:1 ends here

// crystal
// - [[https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html][wwPDB Format version 3.3: Crystallographic and Coordinate Transformation Section]]


// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*crystal][crystal:1]]
// Example
// -------
// CRYST1   18.126   18.126    7.567  90.00  90.00 120.00 P6/MMM
// ORIGX1      1.000000  0.000000  0.000000        0.00000
// ORIGX2      0.000000  1.000000  0.000000        0.00000
// ORIGX3      0.000000  0.000000  1.000000        0.00000
// SCALE1      0.055169  0.031852  0.000000        0.00000
// SCALE2      0.000000  0.063704  0.000000        0.00000
// SCALE3      0.000000  0.000000  0.132153        0.00000

// Record Format
// -------------
//  COLUMNS      DATA  TYPE    FIELD          DEFINITION
//  -------------------------------------------------------------
//  1 -  6       Record name   "CRYST1"
//  7 - 15       Real(9.3)     a              a (Angstroms).
//  16 - 24      Real(9.3)     b              b (Angstroms).
//  25 - 33      Real(9.3)     c              c (Angstroms).
//  34 - 40      Real(7.2)     alpha          alpha (degrees).
//  41 - 47      Real(7.2)     beta           beta (degrees).
//  48 - 54      Real(7.2)     gamma          gamma (degrees).
//  56 - 66      LString       sGroup         Space  group.
//  67 - 70      Integer       z              Z value.


named!(get_lattice<&str, Lattice>, do_parse!(
                                      tag!("CRYST1")   >>
    a     : flat_map!(take!(9), sp!(parse_to!(f64)))   >>
    b     : flat_map!(take!(9), sp!(parse_to!(f64)))   >>
    c     : flat_map!(take!(9), sp!(parse_to!(f64)))   >>
    alpha : flat_map!(take!(7), sp!(parse_to!(f64)))   >>
    beta  : flat_map!(take!(7), sp!(parse_to!(f64)))   >>
    gamma : flat_map!(take!(7), sp!(parse_to!(f64)))   >>
            take!(1)                                   >>
    sgroup: take!(11)                                  >>
            read_until_eol                             >>
    (
        {
            let sgroup = sgroup.trim();
            let mut lat = Lattice::from_params(a, b, c, alpha, beta, gamma);

            lat
        }
    )
));

#[test]
fn test_pdb_lattice() {
    let lines = "CRYST1   18.126   18.126    7.567  90.00  90.00 120.00 P6/MMM
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.055169  0.031852  0.000000        0.00000
SCALE2      0.000000  0.063704  0.000000        0.00000
SCALE3      0.000000  0.000000  0.132153        0.00000
ATOM      1  O2  MOL     2      -4.808   4.768   2.469  1.00  0.00           O
ATOM      2  O3  MOL     2      -6.684   6.549   1.983  1.00  0.00           O
ATOM      3 T1   MOL     2      -5.234   6.009   1.536  1.00  0.00          Si1+
";
    let (r, v) = get_lattice(lines).expect("pdb lattice");
}
// crystal:1 ends here

// atom record

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*atom%20record][atom record:1]]
// guess element from columns 55-80
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
fn guess_element<'a>(name: &'a str, remained: Option<&'a str>) -> Option<&'a str> {
    // 1. return element symbol without whitespace
    if let Some(r) = remained {
        if let Some(sym) = r.get(22..24).and_then(|s| Some(s.trim())) {
            if ! sym.is_empty() {
                return Some(sym);
            }
        }
    }

    // 2. check atom name
    // ignore the first char if it is a digit
    if let Some(e) = name.chars().next() {
        if ! e.is_alphabetic() {
            return name.get(1..2);
        }
    }
    return name.get(0..1);
}

#[test]
fn test_guess_element() {
    // case 1: with columns containing element
    let x = guess_element("1CA ", Some("  1.00  0.00      UC1 SI"));
    assert_eq!(Some("SI"), x);
    let x = guess_element("1CA ", Some("  1.00  0.00      UC1  I"));
    assert_eq!(Some("I"), x);

    // case 2: without columns containing element
    let x = guess_element("CA  ", None);
    assert_eq!(Some("C"), x);
    let x = guess_element("1SA  ", None);
    assert_eq!(Some("S"), x);
    let x = guess_element(" N B ", None);
    assert_eq!(Some("N"), x);
    // when the remained is just whitespace
    let x = guess_element(" H   ", Some("                        "));
    assert_eq!(Some("H"), x);
}

// ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474  1.00  0.00      UC1 SI
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
named!(atom_record<&str, (usize, Atom)>, do_parse!(
    // 1-6
               alt!(tag!("ATOM  ") | tag!("HETATM")) >>
    // 7-11
    sn       : flat_map!(take!(5), sp!(parse_to!(usize))) >>
    // 12
               take!(1)                               >>
    // 13-16
    name     : take!(4)                               >>
    // 17
    alt_loc  : take!(1)                               >>
    // 18-20
    res_name : take!(3)                               >>
    // 21
               take!(1)                               >>
    // 22
    chain_id : take!(1)                               >>
    // 23-26
    res_seq  : take!(4)                               >>
    // 27
    icode    : take!(1)                               >>
    // 28-30
               take!(3)                               >>
    // 31-38
    x        : flat_map!(take!(8), sp!(parse_to!(f64)))     >>
    // 39-46
    y        : flat_map!(take!(8), sp!(parse_to!(f64)))     >>
    // 47-54
    z        : flat_map!(take!(8), sp!(parse_to!(f64)))     >>
    remained : opt!(read_until_eol)           >>
    (
        {
            let sym = guess_element(name, remained).unwrap();

            let mut a = Atom::new(sym, [x, y, z]);

            (sn, a)
        }
    )
));

#[test]
fn test_pdb_atom() {
    let line = "ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474\n";
    let line = "HETATM 1632  O1S MID E   5      -6.883   5.767  26.435  1.00 26.56           O \n";
    let x = atom_record(line);
    println!("{:?}", x);
    let (_, (i, a)) = x.unwrap();

    let line2 = format_atom(3, &a);
    println!("{:?}", line);
    println!("{:?}", line2);
}

named!(get_atoms_from<&str, Vec<(usize, Atom)>>, do_parse!(
    atoms: many0!(atom_record)                            >>
    (
        atoms
    )
));

#[test]
fn test_pdb_get_atoms() {
    let lines = "HETATM 1631  S   MID E   5      -5.827   4.782  25.917  1.00 24.57           S
HETATM 1634  C1  MID E   5      -3.761   3.904  27.580  1.00 28.14           C
ATOM   1634  C1  MID E   5      -3.761   3.904  27.580  1.00 28.14           C
HETATM 1641  C8  MID E   5      -2.096   3.018  29.071  1.00 30.82           C\n\n";
    let (_, atoms) = get_atoms_from(lines).expect("pdb atoms");
    assert_eq!(4, atoms.len());
}

fn format_atom(i: usize, a: &Atom) -> String {
    let [x, y, z] = a.position();

    format!(
        "ATOM  {index:>5} {name:<4}{alt_loc:1}{res_name:<3} {chain_id:1}{res_seq:>4}{icode:>1}   {x:-8.3}{y:-8.3}{z:-8.3}  1.00  0.00          {symbol:>2}\n",
        index=i,
        alt_loc=" ",
        res_name="xx",
        name=a.label(),
        chain_id=1,
        res_seq=1,
        icode=" ",
        x = x,
        y = y,
        z = z,
        symbol=a.symbol(),
    )
}
// atom record:1 ends here

// bond record

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*bond%20record][bond record:1]]
use crate::parser::space_token;

// atom index in bond record line
// named!(bond_atom_index<&str, usize>, flat_map!(
//     take!(5),
//     sp!(parse_to!(usize))
// ));

// https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html
//
// COLUMNS       DATA  TYPE      FIELD        DEFINITION
// -------------------------------------------------------------------------
// 1 -  6        Record name    "CONECT"
// 7 - 11        Integer        serial       Atom  serial number
// 12 - 16       Integer        serial       Serial number of bonded atom
// 17 - 21       Integer        serial       Serial  number of bonded atom
// 22 - 26       Integer        serial       Serial number of bonded atom
// 27 - 31       Integer        serial       Serial number of bonded atom
//
// Example
// -------
// CONECT 1179  746 1184 1195 1203
// CONECT 1179 1211 1222
// CONECT 1021  544 1017 1020 1022
// NOTE: Expected to fail if atom index is larger than 9999 since
// neighboring numbers will overlap
named!(bond_record<&str, Vec<(usize, usize)>>, do_parse!(
             tag!("CONECT")                       >>
    current: sp!(unsigned_digit)                  >>
    others : many_m_n!(1, 4, sp!(unsigned_digit)) >>
             sp!(line_ending)                     >>
    (
        {
            let mut pairs = vec![];
            for other in others {
                pairs.push((current, other));
            }

            pairs
        }
    )
));

fn format_bonds(mol: &Molecule) -> String {
    let mut lines = String::new();

    // connectivity
    // FIXME: add new method in molecule
    let mut map = HashMap::new();
    for (i, j, b) in mol.view_bonds() {
        let mut neighbors = map.entry(i).or_insert(vec![]);
        neighbors.push((j, b.order()));
    }
    for (i, a) in mol.view_atoms() {
        if let Some(neighbors) = map.get(&i) {
            let mut line = format!("CONECT{:>5}", i);
            for (j, _) in neighbors {
                line.push_str(&format!("{:>5}", j));
            }
            lines.push_str(&format!("{}\n", line));
        }
    }

    lines
}

#[test]
fn test_pdb_bond_record() {
    let line = "CONECT 1179 1211 1222 \n";
    let (_, x) = bond_record(line)
        .expect("pdb bond record test1");
    assert_eq!(2, x.len());

    let line = "CONECT 2041 2040 2042\n";
    let (_, x) = bond_record(line).unwrap();
    assert_eq!(2, x.len());

    let line = "CONECT 1179  746 11        \n";
    let (r, x) = bond_record(line).unwrap();
    assert_eq!(2, x.len());
}

named!(get_bonds_from<&str, Vec<(usize, usize)>>, do_parse!(
    bonds: many0!(bond_record) >>
    (
        bonds.into_iter().flat_map(|x| x).collect()
    )
));

#[test]
fn test_pdb_get_bonds() {
    let lines = "CONECT 2028 2027 2029
CONECT 2041 2040 2042
CONECT 2043 2042 2044
\n";

    let (_, x) = get_bonds_from(lines)
        .expect("pdb bonds");
    assert_eq!(6, x.len());

    let lines = "CONECT 2028 2027 2029
\n";

    let (_, x) = get_bonds_from(lines)
        .expect("pdb missing bonds");
    assert_eq!(2, x.len());
}
// bond record:1 ends here

// molecule

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*molecule][molecule:1]]
use std::collections::HashMap;

named!(get_molecule<&str, Molecule>, do_parse!(
    // 1. read lattice data
    // locate loattice record
    many_till!(read_until_eol, peek!(
        alt!(
            tag!("CRYST1") |
            alt!(
                tag!("ATOM") |
                tag!("HETATM")
            ))
    )) >>
    lattice: opt!(complete!(get_lattice))                   >>
    // 2. read atoms
    // locate atom records
    many_till!(read_until_eol, peek!(atom_record)) >>
    atoms  : get_atoms_from                                 >>
    // 3. read bonds
    // locate bond records
    many_till!(read_until_eol, alt!(peek!(tag!("CONECT")) |
                                    terminated!(
                                        alt!(
                                            tag!("END") |
                                            tag!("ENDMDL")
                                        ),
                                        line_ending)))      >>
    bonds  : opt!(complete!(get_bonds_from))                >>
    (
        {
            let mut mol = Molecule::new("pdb");
            // table for mapping atom id and atom index
            let mut table = HashMap::new();
            for (i, a) in atoms {
                let n = mol.add_atom(a);
                table.insert(i, n);
            }

            // add bonds
            if let Some(bonds) = bonds {
                let mut bonds = bonds;
                for b in bonds {
                    let n1 = table[&b.0];
                    let n2 = table[&b.1];
                    mol.add_bond(n1, n2, Bond::single());
                }
            }

            mol.lattice = lattice;
            mol
        }
    )
));

#[test]
fn test_pdb_molecule() {
    let lines = "SCALE3      0.000000  0.000000  0.132153        0.00000
ATOM      1  O2  MOL     2      -4.808   4.768   2.469  1.00  0.00           O
ATOM      2  O3  MOL     2      -6.684   6.549   1.983  1.00  0.00           O
ATOM      3 T1   MOL     2      -5.234   6.009   1.536  1.00  0.00          Si1+
ATOM      4  O1  MOL     2      -4.152  10.936   1.688  1.00  0.00           O
ATOM      5  O1  MOL     2      -4.150  10.935   1.688  1.00  0.00           O
ATOM      6  O2  MOL     2      -1.725  11.578   2.469  1.00  0.00           O
ATOM      7  O2  MOL     2      -9.164  10.843   2.469  1.00  0.00           O
ATOM      8 T1   MOL     2      -2.587  10.589   1.536  1.00  0.00          Si1+
ATOM      9 T1   MOL     2      -7.877  10.591   1.536  1.00  0.00          Si1+
ATOM     10  O2  MOL     2      -1.725  -6.548   2.469  1.00  0.00           O
ATOM     11  O3  MOL     2      -2.330  -9.063   1.983  1.00  0.00           O
ATOM     12 T1   MOL     2      -2.587  -7.537   1.536  1.00  0.00          Si1+
ATOM     13  O1  MOL     2      -7.395  -9.064   1.688  1.00  0.00           O
TER     367
CONECT    2    4
CONECT    3    4
END\n";

    let (r, v) = get_molecule(lines).expect("pdb molecule");
    assert_eq!(13, v.natoms());
    assert_eq!(2, v.nbonds());

        let lines = "SCALE3      0.000000  0.000000  0.132153        0.00000
ATOM      1  O2  MOL     2      -4.808   4.768   2.469  1.00  0.00           O
ATOM      2  O3  MOL     2      -6.684   6.549   1.983  1.00  0.00           O
ATOM      3 T1   MOL     2      -5.234   6.009   1.536  1.00  0.00          Si1+
ATOM      4  O1  MOL     2      -4.152  10.936   1.688  1.00  0.00           O
ATOM      5  O1  MOL     2      -4.150  10.935   1.688  1.00  0.00           O
ATOM      6  O2  MOL     2      -1.725  11.578   2.469  1.00  0.00           O
ATOM      7  O2  MOL     2      -9.164  10.843   2.469  1.00  0.00           O
ATOM      8 T1   MOL     2      -2.587  10.589   1.536  1.00  0.00          Si1+
ATOM      9 T1   MOL     2      -7.877  10.591   1.536  1.00  0.00          Si1+
ATOM     10  O2  MOL     2      -1.725  -6.548   2.469  1.00  0.00           O
ATOM     11  O3  MOL     2      -2.330  -9.063   1.983  1.00  0.00           O
ATOM     12 T1   MOL     2      -2.587  -7.537   1.536  1.00  0.00          Si1+
ATOM     13  O1  MOL     2      -7.395  -9.064   1.688  1.00  0.00           O
END\n";

    let (r, v) = get_molecule(lines).expect("pdb molecule no bonds");
    assert_eq!(13, v.natoms());
    assert_eq!(0, v.nbonds());
}

fn format_molecule(mol: &Molecule) -> String {
    if mol.natoms() > 9999 {
        eprintln!("WARNING: PDB format is incapable for large molecule (natoms < 9999)");
    }

    // atoms
    let mut lines = String::from("REMARK Created by gchemol\n");
    for (i, a) in mol.view_atoms() {
        let line = format_atom(i, a);
        lines.push_str(&line);
    }

    // bonds
    if mol.nbonds() > 0 {
        lines.push_str(&format_bonds(&mol));
    }

    lines.push_str("END\n");

    lines
}
// molecule:1 ends here

// chemfile

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*chemfile][chemfile:1]]
pub struct PdbFile();

impl ChemFileLike for PdbFile {
    fn ftype(&self) -> &str {
        "text/pdb"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".pdb", ".ent"]
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        Ok(format_molecule(mol))
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule(chunk)
    }
}

#[test]
fn test_pdb_parse() {
    let file = PdbFile();
    // single molecule with a lattice
    let mols = file.parse("tests/files/pdb/sio2.pdb")
        .expect("parse pdb molecules");
    assert_eq!(1, mols.len());
    assert!(mols[0].lattice.is_some());
}
// chemfile:1 ends here
