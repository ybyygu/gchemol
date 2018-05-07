// [[file:~/Workspace/Programming/gchemol/gchemol.note::ff90f63c-4f42-4a44-8333-59dac76a029f][ff90f63c-4f42-4a44-8333-59dac76a029f]]
use Atom;
use parser::{
    end_of_line,
    unsigned_digit,
    double_s,
    not_space,
};

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
named!(atom_record<&str, Atom>, do_parse!(
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
            Atom::new(e, [x, y, z])
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
    let (r, a) = atom_record(" 3	C3	2.414	0.000	0.000	C.ar	1	BENZENE	0.000").unwrap();
    assert_eq!("C", a.symbol());
}
// ff90f63c-4f42-4a44-8333-59dac76a029f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::e6d75d58-cab8-47f3-85ea-e710192a4a82][e6d75d58-cab8-47f3-85ea-e710192a4a82]]
use parser::{
    alphanumeric,
    space_token,
    take_until_end_of_line,
};

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
named!(get_bonds_from<&str, Vec<(usize, usize, &str)>>, do_parse!(
    tag!("@<TRIPOS>BOND")  >>
        take_until_end_of_line >>
        bonds: many0!(bond_record)    >>
        (bonds)
));

#[test]
fn test_mol2_bonds() {
    let (_, x) = get_bonds_from("@<TRIPOS>BOND
1 1 2 1
2 1 3 1
").unwrap();
    assert_eq!(2, x.len());

    let (_, x) = get_bonds_from("@<TRIPOS>BOND\n").unwrap();
    assert_eq!(0, x.len());
}
// e6d75d58-cab8-47f3-85ea-e710192a4a82 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::13916972-0d19-4b09-807f-e1d45ac3ab2b][13916972-0d19-4b09-807f-e1d45ac3ab2b]]
use Molecule;

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
                  take_until!("@<TRIPOS>")         >>
    bonds       : opt!(get_bonds_from)                 >>
    (
        {
            if natoms != atoms.len() {
                eprintln!("Inconsistency: expected {} atoms, but found {}", natoms, atoms.len());
            }
            let mut mol = Molecule::new(title);
            for a in atoms {
                mol.add_atom(a);
            }

            mol
        }
    )
));

#[test]
fn test_formats_mol2() {
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

    println!("{:?}", x);

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

use nom::IResult;

use errors::*;
use formats::{
    ChemFileLike,
};

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
}

#[test]
fn test_formats_xyz() {
    let file = Mol2File();
    let mols = file.parse("tests/files/mol2/multi-obabel.mol2").unwrap();
    assert_eq!(6, mols.len());
}
// 91746805-686b-489a-b077-2b7182f18e3b ends here
