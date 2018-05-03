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
// e6d75d58-cab8-47f3-85ea-e710192a4a82 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::5239f37d-ccc7-4f8f-9200-ecfff68c5751][5239f37d-ccc7-4f8f-9200-ecfff68c5751]]

// 5239f37d-ccc7-4f8f-9200-ecfff68c5751 ends here
