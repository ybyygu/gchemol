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
