// header

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*header][header:1]]
// parses the following record types in a PDB file:
//
// CRYST
// ATOM or HETATM
// TER
// END or ENDMDL
// CONECT
// header:1 ends here

// base

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*base][base:1]]
use super::*;
// base:1 ends here

// crystal
// # References
// - [[https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html][wwPDB Format version 3.3: Crystallographic and Coordinate Transformation Section]]

// # Example
// CRYST1   18.126   18.126    7.567  90.00  90.00 120.00 P6/MMM
// ORIGX1      1.000000  0.000000  0.000000        0.00000
// ORIGX2      0.000000  1.000000  0.000000        0.00000
// ORIGX3      0.000000  0.000000  1.000000        0.00000
// SCALE1      0.055169  0.031852  0.000000        0.00000
// SCALE2      0.000000  0.063704  0.000000        0.00000
// SCALE3      0.000000  0.000000  0.132153        0.00000

// # Record Format
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

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*crystal][crystal:1]]
named!(read_lattice<&str, Lattice>, do_parse!(
                                      tag!("CRYST1")   >>
    a     : flat_map!(take!(9), sp!(parse_to!(f64)))   >>
    b     : flat_map!(take!(9), sp!(parse_to!(f64)))   >>
    c     : flat_map!(take!(9), sp!(parse_to!(f64)))   >>
    alpha : flat_map!(take!(7), sp!(parse_to!(f64)))   >>
    beta  : flat_map!(take!(7), sp!(parse_to!(f64)))   >>
    gamma : flat_map!(take!(7), sp!(parse_to!(f64)))   >>
            take!(1)                                   >>
    sgroup: take!(11)                                  >>
            read_line                                  >>
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
    let (_, mut v) = read_lattice(lines).expect("pdb lattice");
    let abc = v.lengths();
    assert_eq!(abc[1], 18.126);
}
// crystal:1 ends here

// element
// # guess element from data in columns 55-80
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*element][element:1]]
fn guess_element<'a>(name: &'a str, r: &'a str) -> Option<&'a str> {
    // 1. return element symbol without whitespace
    if let Some(sym) = r.get(22..24).and_then(|s| Some(s.trim())) {
        if ! sym.is_empty() {
            return Some(sym);
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
    // case 1: with columns containing element symbols
    let x = guess_element("1CA ", "  1.00  0.00      UC1 SI");
    assert_eq!(Some("SI"), x);
    let x = guess_element("1CA ", "  1.00  0.00      UC1  I");
    assert_eq!(Some("I"), x);

    // case 2: without columns containing element symbols
    let x = guess_element("CA  ", "");
    assert_eq!(Some("C"), x);
    let x = guess_element("1SA  ", "");
    assert_eq!(Some("S"), x);
    let x = guess_element(" N B ", "");
    assert_eq!(Some("N"), x);
    // when the remained is just whitespace
    let x = guess_element(" H   ", "                        ");
    assert_eq!(Some("H"), x);
}
// element:1 ends here

// atom records
// # Example
// ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474  1.00  0.00      UC1 SI
// # Format
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*atom%20records][atom records:1]]
// Return Atom index (sn) and Atom object
named!(read_atom_record<&str, (usize, Atom)>, do_parse!(
               alt!(tag!("ATOM  ") | tag!("HETATM"))      >> // 1-6
    sn       : flat_map!(take!(5), sp!(parse_to!(usize))) >> // 7-11
               take!(1)                                   >> // 12
    name     : take!(4)                                   >> // 13-16
    alt_loc  : take!(1)                                   >> // 17
    res_name : take!(3)                                   >> // 18-20
               take!(1)                                   >> // 21
    chain_id : take!(1)                                   >> // 22
    res_seq  : take!(4)                                   >> // 23-26
    icode    : take!(1)                                   >> // 27
               take!(3)                                   >> // 28-30
    x        : flat_map!(take!(8), sp!(parse_to!(f64)))   >> // 31-38
    y        : flat_map!(take!(8), sp!(parse_to!(f64)))   >> // 39-46
    z        : flat_map!(take!(8), sp!(parse_to!(f64)))   >> // 47-54
    rest     : read_line                                  >>
    (
        {
            // TODO: take more attributes
            let sym = guess_element(name, rest).unwrap();
            let mut a = Atom::new(sym, [x, y, z]);

            (sn, a)
        }
    )
));

// Render atom in pdb format
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

#[test]
fn test_pdb_atom() {
    let line = "ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474\n";
    let (_, (i, a)) = read_atom_record(line).expect("pdb atom");
    assert_eq!(3, i);
    assert_eq!("S", a.symbol());
    assert_eq!([3.484, 3.484, 3.474], a.position());

    let line = "ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474  1.00  0.00      UC1 SI\n";
    let (_, (i, a)) = read_atom_record(line).expect("pdb atom");
    assert_eq!("Si", a.symbol());

    let line = "HETATM 1632  O1S MID E   5      -6.883   5.767  26.435  1.00 26.56           O \n";
    let (_, (i, a)) = read_atom_record(line).expect("pdb atom");
    assert_eq!(1632, i);
    assert_eq!("O", a.symbol());
    assert_eq!([-6.883, 5.767, 26.435], a.position());

    let line = format_atom(3, &a);
    let (_, (i, b)) = read_atom_record(&line).expect("pdb atom");
    assert_eq!(3, i);
    assert_eq!(a.symbol(), b.symbol());
    assert_eq!(a.position(), b.position());

}

named!(read_atoms<&str, Vec<(usize, Atom)>>, do_parse!(
    atoms: many0!(read_atom_record) >>
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
    let (_, atoms) = read_atoms(lines).expect("pdb atoms");
    assert_eq!(4, atoms.len());
}
// atom records:1 ends here

// bond records
// # References
// - https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html

// # Format
// COLUMNS       DATA  TYPE      FIELD        DEFINITION
// -------------------------------------------------------------------------
//  1 -  6       Record name    "CONECT"
//  7 - 11       Integer        serial       Atom  serial number
// 12 - 16       Integer        serial       Serial number of bonded atom
// 17 - 21       Integer        serial       Serial  number of bonded atom
// 22 - 26       Integer        serial       Serial number of bonded atom
// 27 - 31       Integer        serial       Serial number of bonded atom

// # Example
// CONECT 1179  746 1184 1195 1203
// CONECT 1179 1211 1222
// CONECT 1021  544 1017 1020 1022

// # NOTE
// Expected to fail if atom index is larger than 9999 since neighboring numbers
// will overlap

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*bond%20records][bond records:1]]
named!(read_bond_record<&str, Vec<(usize, usize)>>, sp!(do_parse!(
             tag!("CONECT")                  >>
    current: unsigned_digit                  >>
    others : many1!(unsigned_digit)          >>
             line_ending                     >>
    (
        {
            let mut pairs = vec![];
            for other in others {
                pairs.push((current, other));
            }

            pairs
        }
    )
)));

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
fn test_pdb_read_bond() {
    let line = "CONECT 1179 1211 1222 \n";
    let (_, x) = read_bond_record(line)
        .expect("pdb bond record test1");
    assert_eq!(2, x.len());

    let line = "CONECT 2041 2040 2042\n";
    let (_, x) = read_bond_record(line).unwrap();
    assert_eq!(2, x.len());

    let line = "CONECT 1179  746 11        \n";
    let (r, x) = read_bond_record(line).unwrap();
    assert_eq!(2, x.len());
}

named!(read_bonds<&str, Vec<(usize, usize)>>, do_parse!(
    bonds: many0!(read_bond_record) >>
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

    let (_, x) = read_bonds(lines)
        .expect("pdb bonds");
    assert_eq!(6, x.len());

    let lines = "CONECT 2028 2027 2029
\n";

    let (_, x) = read_bonds(lines)
        .expect("pdb missing bonds");
    assert_eq!(2, x.len());
}
// bond records:1 ends here

// molecule

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*molecule][molecule:1]]
use std::collections::HashMap;

// quick jump to crystal record (optional)
named!(jump1<&str, ()>, do_parse!(
    many_till!(
        read_line,
        peek!(alt!(tag!("CRYST1") | tag!("ATOM") | tag!("HETATM")))) >>
    (
        ()
    )
));

// quick jump to Atom records
named!(jump2<&str, ()>, do_parse!(
    many_till!(
        read_line,
        peek!(alt!(tag!("ATOM") | tag!("HETATM")))) >>
    (
        ()
    )
));

// recognize optional record between Atom and Bond
named!(sep_atom_bond<&str, ()>, do_parse!(
    alt!(tag!("TER") | tag!("END") ) >> read_line >>
    (
        ()
    )
));

// Add atoms/bonds into Molecule `mol`
pub fn parse_atoms_and_bonds<'a>(input: &'a str, mol: &mut Molecule) -> nom::IResult<&'a str, ()> {
    do_parse!(input,
              atoms: read_atoms          >>
                     opt!(sep_atom_bond) >>
              bonds: read_bonds          >>
              (
                  {
                      // bond records lookup table
                      let mut table = HashMap::new();

                      for (i, a) in atoms {
                          let n = mol.add_atom(a);
                          table.insert(i, n);
                      }

                      for b in bonds {
                          let n1 = table[&b.0];
                          let n2 = table[&b.1];
                          mol.add_bond(n1, n2, Bond::single());
                      }
                  }
              )
    )
}

fn read_molecule(input: &str) -> nom::IResult<&str, Molecule> {
    let mut mol = Molecule::new("for pdb file");

    // 1. read lattice record
    let (input, _) = jump1(input)?;
    let (input, line) = peek!(input, read_line)?;
    let input = if line.starts_with("CRYST1") {
        let (input, lat) = read_lattice(input)?;
        mol.set_lattice(lat);

        input
    } else {
        input
    };

    // 2. read atom records and bond records
    let input = if line.starts_with("ATOM") || line.starts_with("HETATM") {
        let (input, _) = parse_atoms_and_bonds(input, &mut mol)?;
        input
    }
    // ignore irrelevant preceeding lines
    else {
        let (input, _) = jump2(input)?;
        let (input, _) = parse_atoms_and_bonds(input, &mut mol)?;
        input
    };

    Ok((input, mol))
}

fn format_molecule(mol: &Molecule) -> String {
    if mol.natoms() > 9999 {
        eprintln!("PDB format is incapable for large molecule (natoms < 9999)");
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

#[test]
fn test_pdb_molecule() {
    let lines = "\
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
END\n\n";

    let (r, v) = read_molecule(lines).expect("pdb molecule");
    assert_eq!(13, v.natoms());
    assert_eq!(2, v.nbonds());

    let lines = "\
REMARK   Created:  2018-10-22T12:36:28Z
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
\n\n\n";

    let (r, v) = read_molecule(&lines).expect("pdb molecule no bonds");
    assert_eq!(13, v.natoms());
    assert_eq!(0, v.nbonds());
    let txt = "\
CRYST1   54.758   54.758   55.584  90.00  90.00  90.00 P 1           1
ATOM      1  SI1 SIO2X   1       1.494   1.494   0.000  1.00  0.00      UC1 SI
ATOM      2  O11 SIO2X   1       1.194   0.514   1.240  1.00  0.00      UC1  O
ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474  1.00  0.00      UC1 SI
ATOM      4  O12 SIO2X   1       3.784   4.464   4.714  1.00  0.00      UC1  O
ATOM      5  SI3 SIO2X   1       0.995   3.983   1.737  1.00  0.00      UC1 SI
ATOM      6  O13 SIO2X   1       1.975   3.683   2.977  1.00  0.00      UC1  O
ATOM      7  SI4 SIO2X   1       3.983   0.995   5.211  1.00  0.00      UC1 SI
ATOM      8  O14 SIO2X   1       3.003   1.295   6.451  1.00  0.00      UC1  O
ATOM      9  O21 SIO2X   1       1.295   3.003   0.497  1.00  0.00      UC1  O
ATOM     10  O22 SIO2X   1       3.683   1.975   3.971  1.00  0.00      UC1  O
ATOM     11  O23 SIO2X   1       0.514   1.194   5.708  1.00  0.00      UC1  O
ATOM     12  O24 SIO2X   1       4.464   3.784   2.234  1.00  0.00      UC1  O
END\n
";
    let (_, v) = read_molecule(&txt).expect("pdb crystal");
    assert_eq!(12, v.natoms());
    assert!(v.lattice.is_some());
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

    fn parse_molecule<'a>(&self, chunk: &'a str) -> nom::IResult<&'a str, Molecule> {
        read_molecule(chunk)
    }
}
// chemfile:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*test][test:1]]
#[test]
fn test_formats_pdb() {
    let pdb = PdbFile();

    // single molecule with a lattice
    let fname = Path::new("tests/files/pdb/sio2.pdb");
    let mols = pdb.parse(fname).expect("parse pdb molecules");
    assert_eq!(1, mols.len());
    assert!(mols[0].lattice.is_some());

    // ASE generated pdb file
    let fname = Path::new("tests/files/pdb/ase.pdb");
    let mols = pdb.parse(fname).expect("parse pdb molecules");
    assert_eq!(2, mols.len());
    assert_eq!(mols[0].natoms(), 16);
    assert_eq!(mols[1].natoms(), 10);

    let filenames = vec![
        // GaussView generated pdb file
        "tests/files/pdb/gview.pdb",
        // Chem3D generated pdb file
        "tests/files/pdb/chem3d.pdb",
        // Discovery Studio generated pdb file
        "tests/files/pdb/ds.pdb",
        // Material Studio generated pdb file
        "tests/files/pdb/ms.pdb",
    ];
    for fname in filenames {
        let mols = pdb.parse(Path::new(fname)).expect("parse pdb molecules");
        assert_eq!(1, mols.len());
        assert_eq!(mols[0].natoms(), 16);
        assert_eq!(mols[0].nbonds(), 14);
    }

    // multiple molecules generated by babel
    let fname = "tests/files/pdb/multi-babel.pdb";
    let mols = pdb.parse(Path::new(fname)).expect("parse pdb molecules");
    assert_eq!(6, mols.len());
    assert_eq!(mols[5].natoms(), 13);

    // zeolite with connectivity
    let fname = "tests/files/pdb/LTL-zeolite-ms.pdb";
    let mols = pdb.parse(Path::new(fname)).expect("parse pdb molecules");
    assert_eq!(1, mols.len());

    let fname = "tests/files/pdb/tyr-33-conf1.pdb";
    let mols = pdb.parse(Path::new(fname)).expect("parse pdb tyr33");
    assert_eq!(1, mols.len());
    assert_eq!(mols[0].natoms(), 103);

    let fname = "tests/files/pdb/1PPC_ligand.pdb";
    let mols = pdb.parse(Path::new(fname)).expect("parse pdb 1ppc");
    assert_eq!(1, mols.len());
    assert_eq!(mols[0].natoms(), 37);
}
// test:1 ends here
