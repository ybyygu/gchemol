// [[file:~/Workspace/Programming/gchemol/gchemol.note::ffdfbdbc-f657-4961-a2d8-0ae6d9b261d8][ffdfbdbc-f657-4961-a2d8-0ae6d9b261d8]]
use parser::{
    end_of_line,
    take_until_end_of_line,
    unsigned_digit,
    space,
};

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
named!(atom_record<&str, ([f64; 3], Option<&str>)>,
    do_parse!(
        // 1-6
              alt!(tag!("ATOM   ") |
                   tag!("HETATM")) >>
        // 7-11
        sn       : take!(5)        >>
        // 12
                   take!(1)        >>
        // 13-16
        name     : take!(4)        >>
        // 17
        alt_loc  : take!(1)        >>
        // 18-20
        res_name : take!(3)        >>
        // 22
        chain_id : take!(1)        >>
        // 23-26
        res_seq  : take!(4)        >>
        // 27
        icode    : take!(1)        >>
        // 28-30
                   take!(3)        >>
        // 31-38
        x        : take!(8)        >>
        // 39-46
        y        : take!(8)        >>
        // 47-54
        z        : take!(8)        >>
        remained : opt!(take_until_end_of_line) >>
        (
            (
                [
                    x.trim().parse().unwrap(),
                    y.trim().parse().unwrap(),
                    z.trim().parse().unwrap(),
                ],
                guess_element(name, remained)
            )
        )
    )
);

#[test]
fn test_pdb_atom() {
    // let line = "ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474  1.00  0.00      UC1 SI\n";
    let line = "ATOM      3  SI2 SIO2X   1       3.484   3.484   3.474\n";
    let x = atom_record(line);
    println!("{:?}", x);
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
         take_until_end_of_line               >>
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

#[test]
fn test_pdb_bond_record() {
    let line = "CONECT 1179 1211 1222         ";
    let (_, x) = bond_record(line).unwrap();
    assert_eq!(2, x.len());

    let line = "CONECT 2041 2040 2042";
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
CONECT 2043 2042 2044            ";
    let (_, x) = bonds(lines).unwrap();
    assert_eq!(6, x.len());
}
// 0394bb2f-e054-4466-a090-04a0fbe69e03 ends here
