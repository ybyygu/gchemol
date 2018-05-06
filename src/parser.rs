// [[file:~/Workspace/Programming/gchemol/gchemol.note::85054519-f2d5-4c63-994a-78bbe4f9a30f][85054519-f2d5-4c63-994a-78bbe4f9a30f]]
#![macro_use]
use std::fmt::Debug;
use nom::IResult;
pub use nom::{
    alphanumeric,
    eol,
    double_s,
    is_digit,
    digit,
    line_ending,
    space,
    not_line_ending,
    alpha,
};

/// whitespace including one or more spaces or tabs
named!(pub space_token<&str, &str>, eat_separator!(&b" \t"[..]));
macro_rules! sp (
    ($i:expr, $($args:tt)*) => (
        {
            sep!($i, space_token, $($args)*)
        }
    )
);

/// not any whitespace character
/// will not consume "\n" character
named!(pub not_space<&str, &str>, is_not!(" \t\r\n"));

fn dump<T: Debug>(res: IResult<&str,T>) {
    match res {
        IResult::Done(rest, value) => {println!("Done {:?} {:?}",rest,value)},
        IResult::Error(err) => {println!("Err {:?}",err)},
        IResult::Incomplete(needed) => {println!("Needed {:?}",needed)}
    }
}

pub fn end_of_line(input: &str) -> IResult<&str, &str> {
    alt!(input, eof!() | line_ending)
}

// -1, 0, 1, 2, ...
named!(pub signed_digit<&str, isize>,
       map_res!(
           recognize!
               (
                   pair!(opt!(alt!(char!('+') | char!('-'))), digit)
               ),
           str::parse
       )
);

#[test]
fn test_nom_signed_digit() {
    let x = signed_digit("-1");
    let x = signed_digit("1");
}

named!(pub unsigned_digit<&str, usize>,
    map_res!(
        digit,
        str::parse
    )
);
// 85054519-f2d5-4c63-994a-78bbe4f9a30f ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c755a2a3-458f-42b1-aeb4-7c89071491ef][c755a2a3-458f-42b1-aeb4-7c89071491ef]]
named!(pub take_until_end_of_line<&str, &str>,
    terminated!(
        not_line_ending,
        end_of_line
    )
);

#[test]
fn test_nom_one_line() {
    let x = take_until_end_of_line("this is the end\nok\n").unwrap();
}

named!(pub digit_one_line<&str, usize>,
    map_res!(
        terminated!(digit, end_of_line),
        str::parse
    )
);

#[test]
fn test_nom_digit_one_line() {
    let (_, n) = digit_one_line("123\n").unwrap();
    assert_eq!(123, n);
    let (_, n) = digit_one_line("123").unwrap();
    assert_eq!(123, n);
    let (_, n) = digit_one_line("123\nnext line").unwrap();
    assert_eq!(123, n);
}
// c755a2a3-458f-42b1-aeb4-7c89071491ef ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::f1450c5e-a4b0-4dff-a548-af7fe8e52660][f1450c5e-a4b0-4dff-a548-af7fe8e52660]]
use Atom;
use Molecule;

/// Consume three float numbers separated by one or more spaces
/// Return position array
named!(pub xyz_array<&str, [f64; 3]>,
       do_parse!
       (
        x: double_s >>
        space       >>
        y: double_s >>
        space       >>
        z: double_s >>
        ([x, y, z])
       )
);

#[test]
fn test_nom_xyz_array() {
    let (_, x) = xyz_array("-11.4286  1.7645  0.0000").unwrap();
    assert_eq!(x[2], 0.0);

    let (_, x) = xyz_array("-11.4286  1.7645  0.0000\n").unwrap();
    assert_eq!(x[2], 0.0);

    let (_, x) = xyz_array("-11.4286\t1.7E-5  0.0000").unwrap();
    assert_eq!(x[2], 0.0);
}

named!(symbol_and_position<&str, (&str, [f64; 3])>,
    do_parse!(
        symbol   : alphanumeric    >>
                   space           >>
        position : xyz_array       >>
        (symbol, position)
    )
);

#[test]
fn test_nom_symbol_and_position() {
    let (_, (sym, position)) = symbol_and_position("C1 -11.4286 1.7645  0.0000\n").unwrap();
    assert_eq!(sym, "C1");

    let (_, (sym, position)) = symbol_and_position("1\t -11.4286 1.7645  0.0000\n").unwrap();
    assert_eq!(sym, "1");
}

named!(xyz_atom<&str, Atom>,
    do_parse!(
               opt!(space)         >>
        parts: symbol_and_position >>
               opt!(space)         >>
               end_of_line         >>
        (
            Atom::new(parts.0, parts.1)
        )
    )
);

#[test]
fn test_nom_xyz_atom() {
    let (_, x) = xyz_atom("C -11.4286 -1.3155  0.0000 ").unwrap();
    assert_eq!("C", x.symbol());
}

fn new_molecule(title: &str, atoms: Vec<Atom>) -> Molecule {
    let mut mol = Molecule::new(title);
    for a in atoms {
        mol.add_atom(a);
    }

    mol
}

named!(xyz_molecule<&str, Molecule>,
    do_parse!(
                opt!(space)                >>
        n     : digit_one_line             >>
                opt!(space)                >>
        title : take_until_end_of_line                   >>
        atoms : many0!(xyz_atom) >>
        (
            new_molecule(title, atoms)
        )
    )
);

named!(xyz_molecule_many<&str, Vec<Molecule>>,
       do_parse!(
           mols : many0!(xyz_molecule) >>
           (mols)
       )
);

#[test]
fn test_nom_xyz_molecule() {
    let txt = "12
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
12
Molecule 2  0.000000
C -4.9186  2.2028  0.0000
C -3.5849  1.4328  0.0000
C -3.5849 -0.1072  0.0000
C -4.9186 -0.8772  0.0000
C -6.2523 -0.1072  0.0000
C -6.2523  1.4328  0.0000
H -4.9186  3.2928  0.0000
H -2.6410  1.9778  0.0000
H -2.6410 -0.6522  0.0000
H -4.9186 -1.9672  0.0000
H -7.1963 -0.6522  0.0000
H -7.1963  1.9778  0.00003
C -11.4286  1.7645  0.0000
C -10.0949  0.9945  0.0000
C -10.0949 -0.5455  0.0000
";
    let (_, mols) = xyz_molecule_many(txt).unwrap();
    assert_eq!(2, mols.len());
}

named!(maybe_f64_s<&str, f64>,
       map_res!(is_not_s!(" \t\n"), str::parse)
);
// f1450c5e-a4b0-4dff-a548-af7fe8e52660 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c41ceaa0-01c0-4848-b1ea-3f77e0a3e0fc][c41ceaa0-01c0-4848-b1ea-3f77e0a3e0fc]]
named!(gjf_link0_one_line<&str, (&str)>,
   do_parse!
       (opt!(space)          >>
        char!('%')           >>
        cmd: not_line_ending >>
        line_ending          >>
        (cmd)
       )
);

#[test]
fn test_nom_gjf_link0_one_line() {
    let (_, cmd) = gjf_link0_one_line("%Mem=64MB\n").unwrap();
    assert_eq!("Mem=64MB", cmd);

    let (_, cmd) = gjf_link0_one_line(" %save\n").unwrap();
    assert_eq!("save", cmd);
}

named!(gjf_link0_section<&str, Vec<&str>>,
       many0!(gjf_link0_one_line)
);

#[test]
fn test_nom_gjf_link0_section() {
    let x = gjf_link0_section("%chk=C5H12.chk
%nproc=8
%mem=5GB
#p opt freq=noraman nosymm B3LYP/6-31+G** test geom=connect
");
}

named!(gjf_route_section<&str, &str>,
   do_parse!
       (opt!(space)          >>
        char!('#')           >>
        cmd: ws!(take_until!("\n\n")) >>
        (cmd)
       )
);

#[test]
fn test_nom_gjf_route_section() {
    let x = gjf_route_section("#opt freq=noraman nosymm B3LYP/6-31+G** test geom=connect\n\n");
    let x = gjf_route_section("#p opt freq=noraman nosymm\n B3LYP/6-31+G** test geom=connect\n\n");
}

named!(gjf_title_section<&str, &str>,
   ws!(take_until!("\n\n"))
);
// c41ceaa0-01c0-4848-b1ea-3f77e0a3e0fc ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::dac5abf9-43a6-40a7-bf33-8338e106f738][dac5abf9-43a6-40a7-bf33-8338e106f738]]
use std::collections::HashMap;
use std::iter::FromIterator;

#[derive(Debug)]
struct GaussianAtom<'a> {
    element_label : &'a str,
    mm_type       : Option<&'a str>,
    mm_charge     : Option<f64>,
    frozen_code   : Option<isize>,
    position      : [f64; 3],
    properties    : HashMap<&'a str, &'a str>,
    oniom_layer   : Option<&'a str>,
}

impl<'a> Default for GaussianAtom<'a> {
    fn default() -> Self {
        GaussianAtom {
            element_label : "C",
            mm_type       : None,
            mm_charge     : None,
            frozen_code   : None,
            properties    : HashMap::new(),
            position      : [0.0; 3],
            oniom_layer   : None,
        }
    }
}

// fragment=1
named!(gjf_atom_property<&str, (&str, &str)>,
   do_parse!
       (
           param: alphanumeric >>
           char!('=') >>
           value: alphanumeric >>
           (
               (param, value)
           )
       )
);

#[test]
fn test_nom_gjf_atom_property() {
    let (_, (param, value)) = gjf_atom_property("fragment=1").unwrap();
    assert_eq!("fragment", param);
    assert_eq!("1", value);
}

// (fragment=1,iso=13,spin=3)
named!(gjf_atom_properties<&str, HashMap<&str, &str>>,
   do_parse!
       (
           tag!("(") >>
           properties: separated_list!(tag!(","), gjf_atom_property) >>
           tag!(")") >>
           (
               HashMap::from_iter(properties.into_iter())
           )
       )
);

// MM parameters, such as atom type and partial charge
// -CA--0.25
named!(gjf_atom_mm_params<&str, (&str, Option<f64>)>,
   do_parse!
       (
           tag!("-") >>
           mm_type: alphanumeric >>
           mm_charge: opt!(preceded!(tag!("-"), double_s)) >>
           (
               (mm_type, mm_charge)
           )
       )
);

#[test]
fn test_nom_gjf_atom_mm_params() {
    let (_, (mm_type, mm_charge)) = gjf_atom_mm_params("-CA--0.25").unwrap();
    assert_eq!("CA", mm_type);
    assert_eq!(Some(-0.25), mm_charge);

    let (_, (mm_type, mm_charge)) = gjf_atom_mm_params("-CA ").unwrap();
    assert_eq!("CA", mm_type);
    assert_eq!(None, mm_charge);
}

// C-CA--0.25   0   -4.703834   -1.841116   -0.779093 L
// C-CA--0.25   0   -3.331033   -1.841116   -0.779093 L H-HA-0.1  3
named!(gjf_atom_oniom_params<&str, (&str, Option<&str>)>,
       do_parse!
       (
           opt!(space)                  >>
           layer: alpha                  >>
           link: opt!(take_until!("\n")) >>
           (
               (layer, link)
           )
       )
);


#[test]
fn test_nom_gjf_atom_oniom_params() {
    let x = gjf_atom_oniom_params("L\n").unwrap();
    let (_, (layer, _)) = gjf_atom_oniom_params("L H-Ha-0.1 3\n").unwrap();
    assert_eq!("L", layer);
}

// How about this: C-CA--0.25(fragment=1,iso=13,spin=3) 0 0.0 1.2 3.4 H H-H_
named!(gjf_atom_line<&str, GaussianAtom>,
   do_parse!
       (
                           opt!(space)                 >>
           element_label : alphanumeric                >>
           mm_params     : opt!(gjf_atom_mm_params)    >>
           properties    : opt!(gjf_atom_properties)   >>
                           space                       >>
           frozen        : opt!(terminated!(
                                signed_digit,
                                space
                                ))                     >>
           position      : xyz_array                   >>
           oniom         : opt!(gjf_atom_oniom_params) >>
                           opt!(space)                 >>
                           end_of_line                 >>
           (
               GaussianAtom {
                   element_label: element_label,
                   mm_type: mm_params.and_then(|x| Some(x.0)),
                   mm_charge: mm_params.and_then(|x| x.1),
                   frozen_code: frozen,
                   position: position,
                   properties: if properties.is_some() {properties.unwrap()} else {HashMap::new()},
                   oniom_layer: oniom.and_then(|x| Some(x.0)),
                   ..Default::default()
               }
           )
       )
);

#[test]
fn test_gjf_atom_line() {
    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3)  0.00   0.00   0.00 L H-HA-0.1  3\n\n").unwrap();
    assert_eq!(ga.oniom_layer, Some("L"));

    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3) -1 0.00   0.00   0.00 L \n\n").unwrap();
    assert_eq!(ga.oniom_layer, Some("L"));
    assert_eq!(ga.frozen_code, Some(-1));

    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3) -1 0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.frozen_code, Some(-1));

    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3)  0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.mm_charge, Some(-0.25));
    assert_eq!(ga.mm_type, Some("CA"));
    assert_eq!(ga.properties["fragment"], "1");

    let (_, ga) = gjf_atom_line("C-CA--0.25  0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.mm_charge, Some(-0.25));
    assert_eq!(ga.mm_type, Some("CA"));

    let (_, ga) = gjf_atom_line("C-CA  0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.mm_type, Some("CA"));
    assert_eq!(ga.mm_charge, None);

    let (_, ga) = gjf_atom_line(" C12(fragment=1)  0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.properties["fragment"], "1");
    assert_eq!(ga.position[0], 0.0);
}
// dac5abf9-43a6-40a7-bf33-8338e106f738 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::d03ec7e2-6cc0-475f-8fbc-d140db9ee4b2][d03ec7e2-6cc0-475f-8fbc-d140db9ee4b2]]
// Connectivity lines like this:
// 1 2 1.0 3 1.0 4 1.0 5 1.0
//     2
//     3

named!(gjf_bond_pair<&str, (&str, f64)>,
    do_parse!(
            space    >>
        n:  digit    >>
            space    >>
        o:  double_s >>
        (n, o)
    )
);

fn build_bonds<'a>(index1: &'a str, others: Vec<(&'a str, f64)>) -> Vec<(&'a str, &'a str, f64)> {
    let mut bonds = vec![];
    for (index2, order) in others {
        bonds.push((index1, index2, order));
    }

    bonds
}

named!(gjf_connect_line<&str, Vec<(&str, &str, f64)>>,
    do_parse!
    (
                 opt!(space)           >>
        n:       digit                 >>
        others:  many0!(gjf_bond_pair) >>
                 opt!(space)           >>
                 end_of_line           >>
        (
            build_bonds(n, others)
        )
    )
);

#[test]
fn test_nom_gjf_connectivity() {
    let (_, x) = gjf_connect_line(" 1 2 1.0 3 1.0 4 1.0 5 1.0\n").unwrap();
    assert_eq!(4, x.len());
}
// d03ec7e2-6cc0-475f-8fbc-d140db9ee4b2 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::2357368c-ab7e-4eb3-96a8-a8e4aba19bac][2357368c-ab7e-4eb3-96a8-a8e4aba19bac]]
named!(gjf_molecule<&str, &str>,
    do_parse!
    (
        gjf_link0_section >>
        gjf_route_section >>
        title: gjf_title_section >>
        take_until_end_of_line >>
        ("a")
    )
);
// 2357368c-ab7e-4eb3-96a8-a8e4aba19bac ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ac94be65-37d4-4ea4-b839-e1875c9d45e1][ac94be65-37d4-4ea4-b839-e1875c9d45e1]]
#[test]
fn test_nom_gaussian_input() {
    let txt = "%chk=C5H12.chk
%nproc=8
%mem=5GB
#p opt freq=noraman nosymm B3LYP/6-31+G** test geom=connect

Title Card Required

0 1
 C(Fragment=1)    0         -1.29639700         -0.54790000         -0.04565800 L
 H(Fragment=1)    0         -0.94903500         -1.58509500         -0.09306500 L
 H(Fragment=1)    0         -0.93491200         -0.03582400         -0.94211800 L
 H(Fragment=1)    0         -2.39017900         -0.56423400         -0.09090100 L
 C(Fragment=2)    0         -0.80594100          0.14211700          1.23074100 L
 H(Fragment=2)    0         -1.13863100          1.18750400          1.23095600 L
 H(Fragment=2)    0          0.29065200          0.17299300          1.23044100 L
 C(Fragment=3)    0         -1.29047900         -0.54051400          2.51480900 L
 H(Fragment=3)    0         -0.95681000         -1.58684300          2.51564600 L
 H(Fragment=3)    0         -2.38822700         -0.57344700          2.51397600 L
 C(Fragment=4)    0         -0.80793200          0.14352700          3.79887800 L
 H(Fragment=4)    0         -1.14052400          1.18894500          3.79694800 L
 H(Fragment=4)    0          0.28866200          0.17430200          3.80089400 L
 C(Fragment=5)    0         -1.30048800         -0.54500600          5.07527000 L
 H(Fragment=5)    0         -0.95334000         -1.58219500          5.12438500 L
 H(Fragment=5)    0         -0.94034400         -0.03198400          5.97172800 L
 H(Fragment=5)    0         -2.39434000         -0.56114800          5.11880700 L

 1 2 1.0 3 1.0 4 1.0 5 1.0
 2
 3
 4
 5 6 1.0 7 1.0 8 1.0
 6
 7
 8 9 1.0 10 1.0 11 1.0
 9
 10
 11 12 1.0 13 1.0 14 1.0
 12
 13
 14 15 1.0 16 1.0 17 1.0
 15
 16
 17

--Link1--
%chk=C5H12.chk
#p geom=chk
";

    let x = gjf_molecule(txt);
}
// ac94be65-37d4-4ea4-b839-e1875c9d45e1 ends here
