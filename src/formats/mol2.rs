// [[file:~/Workspace/Programming/gchemol/gchemol.note::ff90f63c-4f42-4a44-8333-59dac76a029f][ff90f63c-4f42-4a44-8333-59dac76a029f]]
use super::*;

// #[derive(Debug)]
// struct AtomRecord<'a> {
//     id: usize,
//     name: &'a str,
//     position: [f64; 3],
//     mm_type: &'a str,
//     subst_id: Option<usize>,
//     subst_name: Option<&'a str>,
//     charge: Option<f64>,
// }

// Sample record
// -------------
// 3	C3	2.414	0.000	0.000	C.ar	1	BENZENE	0.000
//
// Format
// ------
// atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
//
named!(atom_record<&str, (usize, Atom)>, do_parse!(
    // atom index
    id        : sp!(unsigned_digit)                    >>
    // atom name
    name      : sp!(not_space)                         >>
    // cartesian coordinates
    x         : sp!(double_s)                          >>
    y         : sp!(double_s)                          >>
    z         : sp!(double_s)                          >>
    mm_type   : sp!(not_space)                         >>
    // substructure and partial charge, which could be omitted
    optional  : opt!(complete!(atom_subst_and_charge)) >>
    // ignore the others
                take_until_end_of_line                 >>
    (
        {
            let e = mm_type.split(|x| x == '.').next().unwrap();
            let mut a = Atom::new(e, [x, y, z]);
            a.set_label(name.trim());
            // a.properties.set("mm_type", mm_type);
            (id, a)
        }
    )
));

named!(atom_subst_and_charge<&str, (usize,
                                    &str,
                                    Option<f64>)>, do_parse!(
    subst_id   : sp!(unsigned_digit)            >>
        subst_name : sp!(not_space)                 >>
        charge     : opt!(complete!(sp!(double_s))) >>
        (
            {
                (
                    subst_id,
                    subst_name,
                    charge
                )
            }
        )
));

#[test]
fn test_mol2_atom() {
    let (r, (_, a)) = atom_record(" 3	C3	2.414	0.000	0.000	C.ar	1	BENZENE	0.000").unwrap();
    assert_eq!("C", a.symbol());
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

/// simple translation without considering the bonding pattern
/// http://www.sdsc.edu/CCMS/Packages/cambridge/pluto/atom_types.html
/// I just want material studio happy to accept my .mol2 file
fn get_atom_type(atom: &Atom) -> &str {
    match atom.symbol() {
        "C" => "C.3",
        "P" => "P.3",
        "Co" => "Co.oh",
        "Ru" => "Ru.oh",
        "O" => "O.2",
        "N" => "N.3",
        "S" => "S.2",
        "Ti" => "Ti.oh",
        "Cr" => "Cr.oh",
        _ => atom.symbol(),
    }
}
// ff90f63c-4f42-4a44-8333-59dac76a029f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::e6d75d58-cab8-47f3-85ea-e710192a4a82][e6d75d58-cab8-47f3-85ea-e710192a4a82]]
// Sample record
// -----------
// 12	6	12	1
// 	6	5	6	ar
// 5	4	9	am	BACKBONE
//
// Format
// ------
// bond_id origin_atom_id target_atom_id bond_type [status_bits]
//
named!(bond_record<&str, (usize, usize, &str)>, do_parse!(
        sp!(unsigned_digit)    >>
    n1: sp!(unsigned_digit)    >>
    n2: sp!(unsigned_digit)    >>
    bo: sp!(alphanumeric)             >>
        take_until_end_of_line >>
    (
        {
            (n1, n2, bo)
        }
    )
));


#[test]
fn test_mol2_bond() {
    let x = bond_record("1	1	2	ar ");
    let x = bond_record("1	1	2	1 BACKBONE\n");
}

// Sample
// ------
// @<TRIPOS>BOND
// 1 1 2 1
// 2 1 3 1
named!(get_bonds_from<&str, Vec<(usize, usize, Bond)>>, do_parse!(
    tag!("@<TRIPOS>BOND")  >>
        take_until_end_of_line >>
        bonds: many0!(bond_record)    >>
        (
            {
                let mut bs = vec![];
                for (i, j, o) in bonds {
                    let bond = match o.to_lowercase().as_ref() {
                        "1"  => Bond::single(),
                        "2"  => Bond::double(),
                        "3"  => Bond::triple(),
                        "ar" => Bond::aromatic(),
                        "am" => Bond::aromatic(),
                        "nc" => Bond::dummy(),
                        "wc" => Bond::partial(), // gaussian view use this
                        _    => Bond::single()
                    };
                    bs.push((i, j, bond));
                }

                bs
            }
        )
));

#[test]
fn test_mol2_bonds() {
    let (_, x) = get_bonds_from("@<TRIPOS>BOND
     1    13    11    1
     2    11    12    1
     3     8     4    1
     4     7     3    1
     5     4     3   ar ").unwrap();
    assert_eq!(5, x.len());

    let (_, x) = get_bonds_from("@<TRIPOS>BOND\n").unwrap();
    assert_eq!(0, x.len());
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
// e6d75d58-cab8-47f3-85ea-e710192a4a82 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::aa5c8cd2-1665-445b-9737-b1c0ab567ffd][aa5c8cd2-1665-445b-9737-b1c0ab567ffd]]
// Format
// ------
// @<TRIPOS>CRYSIN
// cell cell cell cell cell cell space_grp setting
named!(get_lattice_from<&str, Lattice>, do_parse!(
                tag!("@<TRIPOS>CRYSIN") >>
                take_until_end_of_line  >>
    a         : sp!(double_s)           >>
    b         : sp!(double_s)           >>
    c         : sp!(double_s)           >>
    alpha     : sp!(double_s)           >>
    beta      : sp!(double_s)           >>
    gamma     : sp!(double_s)           >>
    space_grp : sp!(unsigned_digit)     >>
    setting   : sp!(unsigned_digit)     >>
                take_until_end_of_line  >>
    (Lattice::from_params(a, b, c, alpha, beta, gamma))
));

#[test]
fn test_mol2_crystal() {
    let (_, x) = get_lattice_from("@<TRIPOS>CRYSIN
12.312000 4.959000 15.876000 90.000000 99.070000 90.000000 4 1").unwrap();
}
// aa5c8cd2-1665-445b-9737-b1c0ab567ffd ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::13916972-0d19-4b09-807f-e1d45ac3ab2b][13916972-0d19-4b09-807f-e1d45ac3ab2b]]
use std::collections::HashMap;

// Format
// ------
// mol_name
// num_atoms [num_bonds [num_subst [num_feat
//                                  [num_sets]]]]
// mol_type
// charge_type
// [status_bits
// [mol_comment]]
//
named!(get_molecule_from<&str, Molecule>, do_parse!(
                  take_until1!("@<TRIPOS>MOLECUL")     >>
                  take_until_end_of_line               >>
    title       : take_until_end_of_line               >>
    natoms      : sp!(unsigned_digit)                  >>
    nbonds      : opt!(sp!(complete!(unsigned_digit))) >>
                  take_until_end_of_line               >>
    mol_type    : take_until_end_of_line               >>
    charge_type : take_until_end_of_line               >>
                  take_until1!("@<TRIPOS>ATOM")        >>
                  take_until_end_of_line               >>
    atoms       : many1!(atom_record)                  >>
                  take_until!("@<TRIPOS>")             >>
    bonds       : opt!(complete!(get_bonds_from))      >>
    lattice     : opt!(complete!(get_lattice_from))    >>
    (
        {
            if natoms != atoms.len() {
                eprintln!("Inconsistency: expected {} atoms, but found {}", natoms, atoms.len());
            }

            let mut mol = Molecule::new(title);

            // assign atoms
            let mut table = HashMap::new();
            for (i, a) in atoms {
                let n = mol.add_atom(a);
                table.insert(i, n);
            }

            // assign bonds, optionally
            if let Some(bonds) = bonds {
                for (i, j, b) in bonds {
                    let ni = table.get(&i).expect(".mol2 file: look up atom in bond record");
                    let nj = table.get(&j).expect(".mol2 file: look up atom in bond record");
                    mol.add_bond(*ni, *nj, b);
                }
            }

            mol.lattice = lattice;
            // println!("nbonds: {:?}", mol.nbonds());

            mol
        }
    )
));

#[test]
fn test_mol2_molecule() {
    let x = get_molecule_from("# created with PyMOL 2.1.0
@<TRIPOS>MOLECULE
Molecule Name
12
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1	  N1	-1.759	-2.546	0.000	N.3	1	UNK0	0.000
2	  H2	-0.759	-2.575	0.000	 H	1	UNK0	0.000
3	  C3	-2.446	-1.270	0.000	C.3	1	UNK0	0.000
4	  H4	-3.071	-1.193	-0.890	 H	1	UNK0	0.000
5	  C5	-3.333	-1.124	1.232	C.3	1	UNK0	0.000
6	  C6	-1.456	-0.114	0.000	C.2	1	UNK0	0.000
7	  H7	-2.720	-1.189	2.131	 H	1	UNK0	0.000
8	  H8	-3.836	-0.157	1.206	 H	1	UNK0	0.000
9	  H9	-4.077	-1.920	1.241	 H	1	UNK0	0.000
10	 O10	-0.219	-0.346	0.000	O.2	1	UNK0	0.000
11	 H11	-2.284	-3.396	0.000	 H	1	UNK0	0.000
12	 H12	-1.811	0.895	0.000	 H	1	UNK0	0.000
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
3 1 11 1
4 3 4 1
5 3 5 1
6 3 6 1
7 5 7 1
8 5 8 1
9 5 9 1
10 6 10 2
11 6 12 1
@<TRIPOS>SUBSTRUCTURE
1	UNK0	1	GROUP	1 ****	UNK");

    let x = get_molecule_from("# created with PyMOL 2.1.0
@<TRIPOS>MOLECULE
Molecule Name
12
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1	  N1	-1.759	-2.546	0.000	N.3	1	UNK0	0.000
2	  H2	-0.759	-2.575	0.000	 H	1	UNK0	0.000
3	  C3	-2.446	-1.270	0.000	C.3	1	UNK0	0.000
4	  H4	-3.071	-1.193	-0.890	 H	1	UNK0	0.000
5	  C5	-3.333	-1.124	1.232	C.3	1	UNK0	0.000
6	  C6	-1.456	-0.114	0.000	C.2	1	UNK0	0.000
7	  H7	-2.720	-1.189	2.131	 H	1	UNK0	0.000
8	  H8	-3.836	-0.157	1.206	 H	1	UNK0	0.000
9	  H9	-4.077	-1.920	1.241	 H	1	UNK0	0.000
10	 O10	-0.219	-0.346	0.000	O.2	1	UNK0	0.000
11	 H11	-2.284	-3.396	0.000	 H	1	UNK0	0.000
12	 H12	-1.811	0.895	0.000	 H	1	UNK0	0.000

@<TRIPOS>SUBSTRUCTURE
1	UNK0	1	GROUP	1 ****	UNK");

}
// 13916972-0d19-4b09-807f-e1d45ac3ab2b ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::91746805-686b-489a-b077-2b7182f18e3b][91746805-686b-489a-b077-2b7182f18e3b]]
use std::str;
use std::fs::File;
use std::io::prelude::*;

/// Tripos Mol2 File Format
///
/// Reference
/// ---------
/// http://tripos.com/tripos_resources/fileroot/pdfs/mol2_format.pdf
/// http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
///
pub struct Mol2File();

impl ChemFileLike for Mol2File {
    fn ftype(&self) -> &str {
        "text/mol2"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".mol2"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule_from(chunk)
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
            let (a, b, c) = lat.lengths();
            let (alpha, beta, gamma) = lat.angles();
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

        Ok(lines)
    }
}

#[test]
fn test_formats_mol2() {
    let file = Mol2File();
    let mols = file.parse("tests/files/mol2/multi-obabel.mol2").unwrap();
    assert_eq!(6, mols.len());
    let mols = file.parse("tests/files/mol2/LTL-crysin-ds.mol2").unwrap();
    assert_eq!(1, mols.len());
    assert!(mols[0].lattice.is_some());
}
// 91746805-686b-489a-b077-2b7182f18e3b ends here
