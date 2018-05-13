// [[file:~/Workspace/Programming/gchemol/gchemol.note::8f0d2f92-9ab9-43ab-9024-e1c1ff1fe89b][8f0d2f92-9ab9-43ab-9024-e1c1ff1fe89b]]
use super::*;

// MDL SD file format
// Example input
// -------------
//    -1.2940   -0.5496   -0.0457 C   0  0  0  0  0  0  0  0  0  0  0  0
//               x, y, z, symbol = float(line[:10]), float(line[10:20]), float(line[20:30]), line[31:34]
named!(get_atom_from<&str, Atom>, do_parse!(
    x: take!(10) >>
    y: take!(10) >>
    z: take!(10) >>
    s: take!(3) >>
    take_until_end_of_line >>
    (
        {
            let x: f64 = x.trim().parse().unwrap();
            let y: f64 = y.trim().parse().unwrap();
            let z: f64 = z.trim().parse().unwrap();
            let sym = s.trim();

            Atom::new(sym, [x, y, z])
        }
    )
));

#[test]
fn test_sdf_atom() {
    let (_, x) = get_atom_from("   -1.2940   -0.5496   -0.0457 C   0  0  0  0  0  0  0  0  0  0  0  0\n")
        .unwrap();
    assert_eq!("C", x.symbol());
}

//   1  4  1  0  0  0  0
named!(get_bond_from<&str, (&str, &str, Bond)>, do_parse!(
    i: take!(3) >>
    j: take!(3) >>
    b: take!(3) >>
    take_until_end_of_line >>
    (
        {
            let i = i.trim();
            let j = j.trim();
            let b: f64 = b.trim().parse().unwrap();

            (
                i,
                j,
                Bond::new(b),
            )
        }
    )
));

// output atom line in .sdf format
fn format_atom(a: &Atom) -> String {
    let pos = a.position();
    format!(
        "{x:-10.4} {y:-9.4} {z:-9.4} {sym:3} 0  0  0  0  0  0  0  0  0 {index:2}\n",
        x   = pos[0],
        y   = pos[1],
        z   = pos[2],
        sym = a.symbol(),
        index = a.label(),
    )
}

#[test]
fn test_format_sdf_atom() {
    let line = "  -13.5661  206.9157  111.5569 C   0  0  0  0  0  0  0  0  0 12";
    let (_, a) = get_atom_from(line).unwrap();
    let line2 = format_atom(&a);
    assert_eq!(line[..60], line2[..60]);
}

fn format_bond(index1: &str, index2: &str, bond: &Bond) -> String {
    format!(
        "{index1:>3}{index2:>3}{order:3}  0  0  0 \n",
        index1 = index1,
        index2 = index2,
        order = 1
    )
}

#[test]
fn test_sdf_bond() {
    let line = "  6  7  1  0  0  0 ";
    let (_, (index1, index2, bond)) = get_bond_from(line).unwrap();
    let line2 = format_bond(index1, index2, &bond);
    assert_eq!(line[..9], line2[..9]);
}
// 8f0d2f92-9ab9-43ab-9024-e1c1ff1fe89b ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::3c9c404d-6b1c-4a07-a34f-5ab3a2a969ca][3c9c404d-6b1c-4a07-a34f-5ab3a2a969ca]]
fn format_molecule(mol: &Molecule) -> String {
    let mut lines = String::new();

    // molecule name
    lines.push_str(&format!("{}\n", mol.name));
    // software
    lines.push_str("gchemol\n");
    // comment
    lines.push_str("\n");
    // counts line
    let line = format!("{natoms:3}{nbonds:3}  0  0  0  0  0  0  0  0999 V2000 \n",
                       natoms=mol.natoms(),
                       nbonds=mol.nbonds());

    lines.push_str(&line);

    for a in mol.atoms() {
        lines.push_str(&format_atom(a));
    }

    for b in mol.bonds() {
        let (a1, a2)= b.partners(&mol).unwrap();
        let i = a1.label();
        let j = a2.label();
        lines.push_str(&format_bond(&i, &j, &b));
    }

    lines
}
// 3c9c404d-6b1c-4a07-a34f-5ab3a2a969ca ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a81f736e-2244-4cc8-8deb-ce493cbd8325][a81f736e-2244-4cc8-8deb-ce493cbd8325]]
pub struct SdfFile();

impl ChemFileLike for SdfFile {
    fn ftype(&self) -> &str {
        "text/mol"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".sd", ".sdf", ".mol"]
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        Ok(format_molecule(mol))
    }
}
// a81f736e-2244-4cc8-8deb-ce493cbd8325 ends here
