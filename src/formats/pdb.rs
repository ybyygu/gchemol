// [[file:~/Workspace/Programming/gchemol/gchemol.note::fb7c688e-38d6-455b-a8f0-ae54c563f3cf][fb7c688e-38d6-455b-a8f0-ae54c563f3cf]]

// fb7c688e-38d6-455b-a8f0-ae54c563f3cf ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ffdfbdbc-f657-4961-a2d8-0ae6d9b261d8][ffdfbdbc-f657-4961-a2d8-0ae6d9b261d8]]
use super::*;

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

fn format_atom(i: usize, a: &Atom) -> String {
    let [x, y, z] = a.position();

    format!(
        "ATOM  {index:>5} {name:>4}{alt_loc:1}{res_name:3}{chain_id:1}{res_seq:4}{icode:1}   {x:-8.4}{y:-8.4}{z:-8.4}\n",
        index=i,
        alt_loc=1,
        name=a.label(),
        chain_id=1,
        res_seq=1,
        res_name="xx",
        icode=1,
        x = x,
        y = y,
        z = z,
    )
}

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

named!(pdb_atoms<&str, &str>, do_parse!(
    many1!(atom_record) >>
    (
        {
            "a"
        }
    )
));
// ffdfbdbc-f657-4961-a2d8-0ae6d9b261d8 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::0394bb2f-e054-4466-a090-04a0fbe69e03][0394bb2f-e054-4466-a090-04a0fbe69e03]]
#[derive(Debug)]
struct Pair(usize, usize);
use parser::space_token;

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
// Expected to fail if atom index is larger than 9999, which I do not care
named!(bond_record<&str, Vec<Pair>>, do_parse!(
         tag!("CONECT")                       >>
    sns: many_m_n!(2, 5, sp!(unsigned_digit)) >>
         read_until_eol               >>
    (
        {
            let mut pairs = vec![];
            let sn0 = sns[0];
            for &sn in &sns[1..] {
                pairs.push(Pair(sn0, sn));
            }

            pairs
        }
    )
));

fn format_bond(i: usize, j: usize, b: &Bond) -> String {
    unimplemented!()
}

#[test]
fn test_pdb_bond_record() {
    let line = "CONECT 1179 1211 1222         \n";
    let (_, x) = bond_record(line).unwrap();
    assert_eq!(2, x.len());

    let line = "CONECT 2041 2040 2042\n";
    let (_, x) = bond_record(line).unwrap();
    assert_eq!(2, x.len());

    let line = "CONECT 1179  746 11        \n";
    let (r, x) = bond_record(line).unwrap();
    assert_eq!(2, x.len());
}

named!(bonds<&str, Vec<Pair>>, do_parse!(
    bonds: many0!(bond_record) >>
    (
        bonds.into_iter().flat_map(|x| x).collect()
    )
));

#[test]
fn test_pdb_bonds() {
    let lines = "CONECT 2028 2027 2029
CONECT 2041 2040 2042
CONECT 2043 2042 2044 \n\n";
    let (_, x) = bonds(lines)
        .expect("pdb bonds");
    assert_eq!(6, x.len());
}
// 0394bb2f-e054-4466-a090-04a0fbe69e03 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b5bad50c-2b36-4d09-9623-d487f9e2333b][b5bad50c-2b36-4d09-9623-d487f9e2333b]]
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
            many_till!(read_until_eol, tag!("CRYST1")) >>
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
    let lines = "REMARK   Materials Studio PDB file
CRYST1   18.126   18.126    7.567  90.00  90.00 120.00 P6/MMM
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
// b5bad50c-2b36-4d09-9623-d487f9e2333b ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::10e26e11-ab6c-4e7d-884a-6a0e98c8d08f][10e26e11-ab6c-4e7d-884a-6a0e98c8d08f]]
use std::collections::HashMap;

named!(get_molecule_from<&str, Molecule>, do_parse!(
    pair : many_till!(read_until_eol, atom_record) >>
    atoms: many0!(atom_record)                     >>
    bpair: many_till!(read_until_eol, bond_record) >>
    bonds: many0!(bond_record)                     >>
    (
        {
            // insert the first atom record
            let mut atoms = atoms;
            let (_, ar) = pair;
            atoms.insert(0, ar);

            let mut mol = Molecule::new("pdb");
            // table for mapping atom id and atom index
            let mut table = HashMap::new();
            for (i, a) in atoms {
                let n = mol.add_atom(a);
                table.insert(i, n);
            }

            // add bonds
            let mut bonds = bonds;
            let (_, br) = bpair;
            bonds.insert(0, br);
            for bs in bonds {
                for b in bs {
                    let n1 = table[&b.0];
                    let n2 = table[&b.1];
                    mol.add_bond(n1, n2, Bond::single());
                }
            }

            mol
        }
    )
));

#[test]
fn test_pdb_molecule() {
    let lines = "SCALE2      0.000000  0.063704  0.000000        0.00000
SCALE3      0.000000  0.000000  0.132153        0.00000
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
    let (r, v) = get_molecule_from(lines).expect("pdb molecule");
    println!("{:?}", (r, v));
}

fn format_molecule(mol: &Molecule) -> String {
    if mol.natoms() > 9999 {
        eprintln!("WARNING: PDB format is incapable for large molecule (natoms < 9999)");
    }

    let mut lines = String::from("REMARK Created by gchemol\n");
    for (i, a) in mol.view_atoms() {
        let line = format_atom(i, a);
        lines.push_str(&line);
    }

    for (i, j, b) in mol.view_bonds() {
        let line = format_bond(i, j, b);
        lines.push_str(&line);
    }

    lines.push_str("END\n");

    lines
}
// 10e26e11-ab6c-4e7d-884a-6a0e98c8d08f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::f026618f-fdc0-4590-a2b0-211078c30e14][f026618f-fdc0-4590-a2b0-211078c30e14]]
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
        get_molecule_from(chunk)
    }
}
// f026618f-fdc0-4590-a2b0-211078c30e14 ends here
