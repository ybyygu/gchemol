// base

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*base][base:1]]
/// Tripos Mol2 File Format
///
/// Reference
/// ---------
/// http://tripos.com/tripos_resources/fileroot/pdfs/mol2_format.pdf
/// http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
///

use super::*;
// base:1 ends here

// atoms
// # Sample record
// @<TRIPOS>ATOM
//       1 O1            0.000906    8.302448    1.688198 O.3      1 SUBUNIT   -0.0000
//       2 O2           -1.779973    6.533331    2.469112 O.3      1 SUBUNIT    0.0000
//       3 O3           -2.514076    9.013548    1.982554 O.3      1 SUBUNIT   -0.0000
//       4 O4           -1.818038    7.434372    0.000000 O.3      1 SUBUNIT   -0.0000
//       5 O5           -2.534921    4.390612    3.783500 O.3      1 SUBUNIT   -0.0000
//       6 O6            0.000000    5.111131    3.783500 O.3      1 SUBUNIT   -0.0000
//       7 T1           -1.528022    7.820533    1.536101 Si       1 SUBUNIT    0.0000
//       8 T2           -1.518959    5.641709    3.783500 Si       1 SUBUNIT    0.0000

// # Format
// atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]


// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*atoms][atoms:1]]
/// Parse Tripos Atom section
named!(read_atoms<&str, Vec<(usize, Atom)>>, preceded!(
    ws!(tag!("@<TRIPOS>ATOM")),
    many1!(read_atom_record)
));

#[test]
fn test_mol2_get_atoms() {
    let lines = "@<TRIPOS>ATOM
      1 N           1.3863   -0.2920    0.0135 N.ar    1  UNL1       -0.2603
      2 N          -1.3863    0.2923    0.0068 N.ar    1  UNL1       -0.2603
      3 C           0.9188    0.9708   -0.0188 C.ar    1  UNL1        0.0456
      4 C          -0.4489    1.2590   -0.0221 C.ar    1  UNL1        0.0456
      5 C          -0.9188   -0.9709    0.0073 C.ar    1  UNL1        0.0456
      6 C           0.4489   -1.2591    0.0106 C.ar    1  UNL1        0.0456
      7 H           1.6611    1.7660   -0.0258 H       1  UNL1        0.0845
      8 H          -0.8071    2.2860   -0.0318 H       1  UNL1        0.0845
      9 H           0.8071   -2.2861    0.0273 H       1  UNL1        0.0845
     10 H          -1.6611   -1.7660    0.0214 H       1  UNL1        0.0845

";
    let (_, atoms) = read_atoms(lines).expect("mol2 atoms");
    assert_eq!(10, atoms.len());
}

named!(read_atom_record<&str, (usize, Atom)>, sp!(do_parse!(
    id        : unsigned_digit         >>     // atom index
    name      : not_space              >>     // atom name
    x         : double                 >>     // cartesian coordinates
    y         : double                 >>
    z         : double                 >>
    emtype    : mm_type                >>     // Element and Atom type
    // substructure and partial charge, which could be omitted
    optional  : opt!(atom_subst_and_charge) >> eol >>
    (
        {
            let (e, mtype) = emtype;
            let mut a = Atom::new(e, [x, y, z]);
            a.set_label(name.trim());
            (id, a)
        }
    )
)));

// parse mol2 atom type. example records:
// C.2, C.3, C.ar
named!(mm_type<&str, (&str, Option<&str>)>, do_parse!(
    symbol: alpha                                    >>
    mtype : opt!(preceded!(tag!("."), alphanumeric)) >>
    (
        (symbol, mtype)
    )
));

#[test]
fn test_mol2_mmtype() {
    let (_, (sym, mtype)) = mm_type("C.ar\n").expect("mol2 atom type");
    assert_eq!("C", sym);
    assert_eq!(Some("ar"), mtype);

    let (_, (sym, mtype)) = mm_type("C.4\n").expect("mol2 atom type 2");
    assert_eq!("C", sym);
    assert_eq!(Some("4"), mtype);

    let (_, (sym, mtype)) = mm_type("C ").expect("mol atom type: missing mm type");
    assert_eq!("C", sym);
    assert_eq!(None, mtype);
}

// substructure id and subtructure name
named!(atom_subst_and_charge<&str, (usize, &str, Option<f64>)>, sp!(do_parse!(
    subst_id   : unsigned_digit            >>
    subst_name : not_space                 >>
    charge     : opt!(double)              >>
    status_bit : opt!(alpha)               >>
    (
        (subst_id, subst_name, charge)
    )
)));

/// simple translation without considering the bonding pattern
/// http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
/// I just want material studio happy to accept my .mol2 file
fn get_atom_type(atom: &Atom) -> &str {
    match atom.symbol() {
        "C"  => "C.3",
        "P"  => "P.3",
        "Co" => "Co.oh",
        "Ru" => "Ru.oh",
        "O"  => "O.2",
        "N"  => "N.3",
        "S"  => "S.2",
        "Ti" => "Ti.oh",
        "Cr" => "Cr.oh",
        _    => atom.symbol(),
    }
}

fn format_atom(a: &Atom) -> String {
    let position = a.position();
    format!("{name:8} {x:-12.5} {y:-12.5} {z:-12.5} {symbol:8} {subst_id:5} {subst_name:8} {charge:-6.4}\n",
            name  = a.label(),
            x = position[0],
            y = position[1],
            z = position[2],
            // FIXME:
            symbol = get_atom_type(a),
            subst_id = 1,
            subst_name = "SUBUNIT",
            charge = 0.0,
    )
}

#[test]
fn test_formats_mol2_atom() {
    let (r, (_, a)) = read_atom_record(" 3	C3	2.414	0.000	0.000	C.ar	1	BENZENE	0.000	DICT\n")
        .expect("mol2 full");
    assert_eq!("C", a.symbol());
    let (r, (_, a)) = read_atom_record(" 3	C3	2.414	0.000	0.000	C.ar	1	BENZENE	0.000\n")
        .expect("mol2 atom: missing status bit");
    assert_eq!("C", a.symbol());
    let (r, (_, a)) = read_atom_record(" 3	C3	2.414	0.000	0.000	C.ar	1	BENZENE \n")
        .expect("mol2 atom: missing partial charge");
    assert_eq!("C", a.symbol());
    let (r, (_, a)) = read_atom_record(" 3	C3	2.414	0.000	0.000	C.ar\n")
        .expect("mol2 atom: missing substructure");
    assert_eq!("C", a.symbol());
}
// atoms:1 ends here

// bond
// # Sample record
// @<TRIPOS>BOND
//   12	6	12	1
//   	6	5	6	ar
//   5	4	9	am	BACKBONE

// # Format
// bond_id origin_atom_id target_atom_id bond_type [status_bits]


// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*bond][bond:1]]
/// Parse Tripos Bond section
named!(read_bonds<&str, Vec<(usize, usize, Bond)>>, do_parse!(
           tag!("@<TRIPOS>BOND")       >> eol >>
    bonds: many0!(read_bond_record)    >>
    (
        bonds
    )
));

#[test]
fn test_mol2_bonds() {
    let lines = "\
@<TRIPOS>BOND
     1    13    11    1
     2    11    12    1
     3     8     4    1
     4     7     3    1
     5     4     3   ar

";

    let (_, x) = read_bonds(lines).expect("mol2 bonds");
    assert_eq!(5, x.len());
}

named!(read_bond_record<&str, (usize, usize, Bond)>, sp!(do_parse!(
        unsigned_digit >>       // bond_id
    n1: unsigned_digit >>       // origin_atom_id
    n2: unsigned_digit >>       // target_atom_id
    bo: alphanumeric   >>       // bond_type
        read_line      >>       // status_bits
    (
        {
            let bond = match bo.to_lowercase().as_ref() {
                "1"  => Bond::single(),
                "2"  => Bond::double(),
                "3"  => Bond::triple(),
                "ar" => Bond::aromatic(),
                "am" => Bond::aromatic(),
                "nc" => Bond::dummy(),
                "wc" => Bond::partial(), // gaussian view use this
                _    => Bond::single()
            };
            (n1, n2, bond)
        }
    )
)));


#[test]
fn test_formats_mol2_bond_record() {
    let (_, (i, j, b)) = read_bond_record("1	1	2	1 BACKBONE\n")
        .expect("mol2 bond: full");
    assert_eq!(BondKind::Single, b.kind);

    let (_, (i, j, b)) = read_bond_record("1	1	2	1\n")
        .expect("mol2 bond: missing status bits");
    assert_eq!(BondKind::Single, b.kind);

    let (_, (i, j, b)) = read_bond_record("1	1	2	ar\n")
        .expect("mol2 bond: aromatic bond type");
    assert_eq!(BondKind::Aromatic, b.kind);
}

fn format_bond_order(bond: &Bond) -> &str {
    match bond.kind {
        BondKind::Single    => "1",
        BondKind::Double    => "2",
        BondKind::Triple    => "3",
        BondKind::Quadruple => "4",
        BondKind::Aromatic  => "ar",
        BondKind::Partial   => "wk", // gaussian view use this
        BondKind::Dummy     => "nc",
    }
}
// bond:1 ends here

// lattice
// # Format
// @<TRIPOS>CRYSIN
// cell cell cell cell cell cell space_grp setting

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*lattice][lattice:1]]
named!(read_lattice<&str, Lattice>, sp!(do_parse!(
                tag!("@<TRIPOS>CRYSIN") >>
                tag!("\n")              >>
    a         : double                  >>
    b         : double                  >>
    c         : double                  >>
    alpha     : double                  >>
    beta      : double                  >>
    gamma     : double                  >>
    space_grp : unsigned_digit          >>
    setting   : unsigned_digit          >>
                read_line               >>
    (Lattice::from_params(a, b, c, alpha, beta, gamma))
)));

#[test]
fn test_formats_mol2_crystal() {
    let txt = "\
@<TRIPOS>CRYSIN
12.312000 4.959000 15.876000 90.000000 99.070000 90.000000 4 1\n";

    let (_, mut x) = read_lattice(txt)
        .expect("mol2 crystal");
    assert_eq!([12.312, 4.959, 15.876], x.lengths());
}
// lattice:1 ends here

// molecule
// # Format
// - mol_name
// - num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
// - mol_type
// - charge_type
// - [status_bits [mol_comment]]

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*molecule][molecule:1]]
fn read_molecule(input: &str) -> IResult<&str, Molecule> {
    let (input, _) = read_lines_until(input, "@<TRIPOS>MOLECULE")?;

    // 1. read meta data
    let (input, (mut mol, counts)) = sp!(input, do_parse!(
                      tag!("@<TRIPOS>MOLECULE") >> eol >>
        title       : read_line                 >>
        counts      : read_usize_many           >>
        mol_type    : read_line                 >>
        charge_type : read_line                 >>
        (
            {
                let mut mol = Molecule::new(title);

                (mol, counts)
            }
        )
    ))?;

    // 2. Assign atoms
    let (input, _) = read_lines_until(input, "@<TRIPOS>ATOM")?;
    let (input, atoms) = read_atoms(input)?;
    let natoms = counts[0];
    if natoms != atoms.len() {
        eprintln!("Inconsistency: expected {} atoms, but found {}", natoms, atoms.len());
    }
    // assign atoms
    let mut table = HashMap::new();
    for (i, a) in atoms {
        let n = mol.add_atom(a);
        table.insert(i, n);
    }

    // 3. Assign bonds (optional)
    let (input, current) = peek_line(input)?;
    let input = if current.starts_with("@<TRIPOS>BOND") {
        let (input, bonds) = read_bonds(input)?;
        for (i, j, b) in bonds {
            let ni = table.get(&i).expect(".mol2 file: look up atom in bond record");
            let nj = table.get(&j).expect(".mol2 file: look up atom in bond record");
            mol.add_bond(*ni, *nj, b);
        }
        input
    } else {
        input
    };

    // 4. Crystal (optional)
    let (input, _) = many_till!(input, read_line, peek!(
        alt!(
            tag!("@<TRIPOS>CRYSIN") |
            tag!("@<TRIPOS>MOLECULE") |
            eof
        )
    ))?;

    let (input, current) = peek_line(input)?;
    let input = if current.starts_with("@<TRIPOS>CRYSIN") {
        let (input, lat) = read_lattice(input)?;
        mol.set_lattice(lat);
        input
    } else {
        input
    };

    Ok((input, mol))
}
// molecule:1 ends here

// chemfile

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*chemfile][chemfile:1]]
pub struct Mol2File();

impl ChemFileLike for Mol2File {
    fn ftype(&self) -> &str {
        "text/mol2"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".mol2"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        read_molecule(chunk)
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        let natoms = mol.natoms();
        let nbonds = mol.nbonds();

        let mut lines = String::new();
        lines.push_str("#	Created by:	gchemol\n");
        lines.push_str("\n");
        lines.push_str("@<TRIPOS>MOLECULE\n");
        lines.push_str(&format!("{}\n", mol.name));

        // atom count, bond numbers, substructure numbers
        lines.push_str(&format!("{:>5} {:>5}\n",
                                natoms,
                                nbonds));
        // molecule type
        lines.push_str("SMALL\n");
        // customed charges
        lines.push_str("USER CHARGES\n");
        // atoms
        lines.push_str("@<TRIPOS>ATOM\n");

        // format atoms
        for (i, a) in mol.view_atoms() {
            let line = format!("{:5} {}", i, format_atom(&a));
            lines.push_str(&line);
        }

        // format bonds
        if nbonds > 0 {
            lines.push_str("@<TRIPOS>BOND\n");
            let mut sn = 1;
            for (i, j, b) in mol.view_bonds() {
                let line = format!("{sn:4} {bond_i:4} {bond_j:4} {bond_type:3}\n",
                                   sn        = sn,
                                   bond_i    = i,
                                   bond_j    = j,
                                   bond_type = format_bond_order(&b));
                lines.push_str(&line);
                sn += 1;
            }
        }

        // format crystal
        if let Some(mut lat) = &mol.lattice {
            lines.push_str("@<TRIPOS>CRYSIN\n");
            let [a, b, c] = lat.lengths();
            let [alpha, beta, gamma] = lat.angles();
            let line = format!("{a:10.4} {b:10.4} {c:10.4} {alpha:5.2} {beta:5.2} {gamma:5.2} {sgrp} 1\n",
                               a     = a,
                               b     = b,
                               c     = c,
                               alpha = alpha,
                               beta  = beta,
                               gamma = gamma,
                               // FIXME: crystal space group
                               sgrp  = 4);

            lines.push_str(&line);
        }

        // final blank line
        lines.push_str("\n");

        Ok(lines)
    }
}
// chemfile:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*test][test:1]]
#[test]
fn test_formats_mol2() {
    let file = Mol2File();

    let mols = file.parse(Path::new("tests/files/mol2/ch3f-dos.mol2")).expect("mol2 ch3f");
    assert_eq!(1, mols.len());

    // when missing final blank line
    // gaussview generated .mol2 file
    let mols = file.parse(Path::new("tests/files/mol2/alanine-gv.mol2")).expect("gv generated mol2 file");
    assert_eq!(1, mols.len());
    let mol = &mols[0];
    assert_eq!(12, mol.natoms());
    assert_eq!(11, mol.nbonds());

    // molecule trajectory
    // openbabel converted .mol2 file
    let mols = file.parse(Path::new("tests/files/mol2/multi-obabel.mol2")).expect("mol2 multi");

    let natoms_expected = vec![16, 10, 16, 16, 16, 13];
    let natoms: Vec<_> = mols.iter().map(|m| m.natoms()).collect();
    assert_eq!(natoms_expected, natoms);

    let nbonds_expected = vec![14, 10, 14, 14, 14, 12];
    let nbonds: Vec<_> = mols.iter().map(|m| m.nbonds()).collect();
    assert_eq!(nbonds_expected, nbonds);
    assert_eq!(6, mols.len());

    // single molecule with a lattice
    // discovery studio generated .mol2 file
    let mols = file.parse(Path::new("tests/files/mol2/LTL-crysin-ds.mol2")).expect("mol2 crysin");
    assert_eq!(1, mols.len());
    assert!(mols[0].lattice.is_some());

}
// test:1 ends here
