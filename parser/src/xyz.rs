// [[file:~/Workspace/Programming/gchemol/parser/combine.note::59a7b68d-4b1f-47fb-b076-85b1d2166407][59a7b68d-4b1f-47fb-b076-85b1d2166407]]
use super::*;

use gchemol::{
    Atom,
    Molecule,
};

/// match xyz coordinates
fn xyz_coordinates<I>() -> impl Parser<Input = I, Output = [f64; 3]>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    (
        signed_float_number(),
        ws(),
        signed_float_number(),
        ws(),
        signed_float_number()
    ).map(|(x, _, y, _, z)| [x, y, z])
}


fn xyz_atom<I>() -> impl Parser<Input = I, Output = Atom>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    (
        optional(ws()).with(many1::<String, _>(char::alpha_num())).message("element type"),
        ws(),
        xyz_coordinates(),
        read_until_eol(),
    ).map(|(sym, _, xyz, _)| {
        Atom::new(sym, xyz)
    })
}

#[test]
fn test_xyz() {
    let line = "-11.4286 \t 1.7645  0.0000";
    let r = xyz_coordinates().easy_parse(line).unwrap();
    assert_eq!(3, r.0.len());

    let (a, _) = xyz_atom().easy_parse("C -11.4286 -1.3155  0.0000\n").unwrap();
    assert_eq!("C", a.symbol());
    let (a, _) = xyz_atom().easy_parse("6 -11.4286 -1.3155  0.0000 \n").unwrap();
    assert_eq!("C", a.symbol());
}

fn xyz_molecule<I>() -> impl Parser<Input = I, Output = Molecule>
where
    I: Stream<Item = char>,
    I::Error: ParseError<I::Item, I::Range, I::Position>,
{
    // (optional(ws()), unsigned_integer(), optional(ws()), eol()).map(|(_, n, _, _)| n).then(|n| {
    //     (
    //         // title
    //         read_until_eol(),
    //         count(n, try(xyz_atom()))
    //     ).map(|(title, atoms)|) {
    //         // 4. construct molecule
    //         let mut mol = Molecule::new(title.trim());
    //         for a in atoms {
    //             mol.add_atom(a);
    //         }

    //         mol
    //     }
    // })

    (
        // line  1: number of atoms
        (optional(ws()), unsigned_integer(), optional(ws()), eol()).map(|(_, n, _, _)| n),
        // line  2: molecule title
        read_until_eol(),
        // lines 3-: atoms
        many1::<Vec<Atom>, _>(try(xyz_atom()))
    ).map(|(n, title, atoms)| {
        // 4. construct molecule
        let mut mol = Molecule::new(title.trim());
        for a in atoms {
            mol.add_atom(a);
        }

        mol
    })
}

#[test]
fn test_xyz_molecule() {
    let txt = "12
  title
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
  ";
    let (mol, _) = xyz_molecule().easy_parse(txt).unwrap();
    assert_eq!(12, mol.natoms());
}
// 59a7b68d-4b1f-47fb-b076-85b1d2166407 ends here
