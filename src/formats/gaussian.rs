// [[file:~/Workspace/Programming/gchemol/gchemol.note::4138ba02-140e-4bdb-8083-74424610b600][4138ba02-140e-4bdb-8083-74424610b600]]
use super::*;
// 4138ba02-140e-4bdb-8083-74424610b600 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c8f2f29e-b23a-4de3-a888-e6cbfee64760][c8f2f29e-b23a-4de3-a888-e6cbfee64760]]
// Gaussian input file
//
// Reference
// ---------
// http://gaussian.com/input/?tabid=0
// c8f2f29e-b23a-4de3-a888-e6cbfee64760 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c41ceaa0-01c0-4848-b1ea-3f77e0a3e0fc][c41ceaa0-01c0-4848-b1ea-3f77e0a3e0fc]]
// sections are separated by a blank line
named!(blank_line<&str, &str>, sp!(line_ending));

named!(gjf_link0<&str, (&str)>, do_parse!(
    opt!(space)              >>
    char!('%')               >>
    cmd: sp!(read_until_eol) >>
    (cmd)
));

#[test]
fn test_link0() {
    let (_, cmd) = gjf_link0("%Mem=64MB\n").unwrap();
    assert_eq!("Mem=64MB", cmd);

    let (_, cmd) = gjf_link0(" %save\n").unwrap();
    assert_eq!("save", cmd);
}

named!(link0_section<&str, Vec<&str>>,
       many0!(gjf_link0)
);

#[test]
fn test_link0_section() {
    let lines = "%chk=C5H12.chk
%nproc=8
%mem=5GB
#p opt freq=noraman nosymm B3LYP/6-31+G** test geom=connect
";

    let (_, link0s) = link0_section(lines).expect("gjf link0 section");
    assert_eq!(3, link0s.len());
}

named!(route_section<&str, String>, do_parse!(
    sp!(char!('#'))                               >>
    parts: many_till!(read_until_eol, blank_line) >>
    (
        {
            let lines = parts.0;
            lines.join(" ")
        }
    )
));

#[test]
fn test_route_section() {
    let lines = "#opt freq=noraman nosymm B3LYP/6-31+G** test geom=connect

";
    let x = route_section(lines).expect("gjf route section");

    let lines = "#p opt freq=noraman nosymm
B3LYP/6-31+G** test geom=connect

";
    let x = route_section(lines).expect("gjf route section multi-lines");
}

named!(title_section<&str, String>, do_parse!(
    parts: many_till!(read_until_eol, blank_line) >>
    (
        {
            let lines = parts.0;
            lines.join(" ")
        }
    )
));
// c41ceaa0-01c0-4848-b1ea-3f77e0a3e0fc ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::dac5abf9-43a6-40a7-bf33-8338e106f738][dac5abf9-43a6-40a7-bf33-8338e106f738]]
use std::collections::HashMap;
use std::iter::FromIterator;

// Specifies the net electric charge (a signed integer) and the spin
// multiplicity (usually a positive integer)
named!(charge_and_spin<&str, Vec<isize>>, separated_list!(
    comma_or_space,
    sp!(signed_digit)
));

named!(charge_and_spin_line<&str, Vec<isize>>, terminated!(
    charge_and_spin,
    sp!(line_ending)
));

#[test]
fn test_charge_and_spin() {
    let (_, x) = charge_and_spin_line("0	1\n").expect("gjf charge & spin");
    assert_eq!(2, x.len());

    let line = " 0 1 , 0 , 1 -3,2 \n";
    let (r, x) = charge_and_spin_line(line).expect("gjf charge & spin");
    assert_eq!(6, x.len());
}

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

// single property entry
// fragment=1
named!(gjf_atom_property<&str, (&str, &str)>, do_parse!(
    param: alphanumeric >>
    char!('=') >>
    value: alphanumeric >>
    (
        {
            (param, value)
        }
    )
));

#[test]
fn test_gjf_atom_property() {
    let (_, (param, value)) = gjf_atom_property("fragment=1 ").unwrap();
    assert_eq!("fragment", param);
    assert_eq!("1", value);
}

// multiple property entries
// (fragment=1,iso=13,spin=3)
named!(gjf_atom_properties<&str, HashMap<&str, &str>>, do_parse!(
    tag!("(")                                           >>
    properties: separated_list!(tag!(","),
                                sp!(gjf_atom_property)) >>
    tag!(")")                                           >>
    (
        HashMap::from_iter(properties.into_iter())
    )
));

#[test]
fn test_gjf_atom_properties() {
    let (_, d) = gjf_atom_properties("(fragment=1,iso=13,spin=3) ")
        .expect("gjf atom properties");
    assert_eq!(3, d.len())
}

fn is_mm_type(chr: char) -> bool {
    is_alphanumeric(chr as u8) || chr == '_'
    // chr != '-' && chr != ' ' && chr != ',' && chr!= '\t' && chr != '\n' && chr != '('
}

// MM parameters, such as atom type and partial charge
// -CA--0.25
named!(gjf_atom_mm_params<&str, (&str, Option<f64>)>, do_parse!(
    tag!("-")                                       >>
    mm_type: take_while!(is_mm_type)                >>
    mm_charge: opt!(preceded!(tag!("-"), double_s)) >>
    (
        (mm_type, mm_charge)
    )
));

#[test]
#[ignore]
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
named!(gjf_atom_line<&str, GaussianAtom>, do_parse!(
                   opt!(space)                 >>
    element_label: alphanumeric                >>
    mm_params    : opt!(gjf_atom_mm_params)    >>
    properties   : opt!(gjf_atom_properties)   >>
                   comma_or_space              >>
    frozen       : opt!(terminated!(
                        signed_digit,
                        comma_or_space
                        ))                     >>
    x            : sp!(double_s)               >>
                   opt!(complete!(tag!(",")))  >>
    y            : sp!(double_s)               >>
                   opt!(complete!(tag!(",")))  >>
    z            : sp!(double_s)               >>
    oniom        : opt!(gjf_atom_oniom_params) >>
                   opt!(comma_or_space)        >>
                   line_ending                 >>
    (
        GaussianAtom {
            element_label: element_label,
            mm_type: mm_params.and_then(|x| Some(x.0)),
            mm_charge: mm_params.and_then(|x| x.1),
            frozen_code: frozen,
            position: [x, y, z],
            properties: if properties.is_some() {properties.unwrap()} else {HashMap::new()},
            oniom_layer: oniom.and_then(|x| Some(x.0)),
            ..Default::default()
        }
    )
));

#[test]
fn test_gjf_atom_line() {
    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3)  0.00   0.00   0.00 L H-HA-0.1  3\n").unwrap();
    assert_eq!(ga.oniom_layer, Some("L"));

    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3) -1 0.00   0.00   0.00 L \n").unwrap();
    assert_eq!(ga.oniom_layer, Some("L"));
    assert_eq!(ga.frozen_code, Some(-1));

    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3) -1 0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.frozen_code, Some(-1));

    let (_, ga) = gjf_atom_line("C-CA--0.25(fragment=1,iso=13,spin=3) ,0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.mm_charge, Some(-0.25));
    assert_eq!(ga.mm_type, Some("CA"));
    assert_eq!(ga.properties["fragment"], "1");

    let (_, ga) = gjf_atom_line("C-CA--0.25  0.00   0.00   0.00\n").unwrap();
    assert_eq!(ga.mm_charge, Some(-0.25));
    assert_eq!(ga.mm_type, Some("CA"));

    let (_, ga) = gjf_atom_line("C-C_3,0.00,0.00,0.00\n").unwrap();
    assert_eq!(ga.mm_type, Some("C_3"));
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
named!(gjf_bond_pair<&str, (usize, f64)>, do_parse!(
        comma_or_space >>
    n:  unsigned_digit >>
        comma_or_space >>
    o:  double_s       >>
    (n, o)
));

// TODO: recover bond order
fn build_bonds(index1: usize, others: Vec<(usize, f64)>) -> Vec<(usize, usize, Bond)> {
    let mut bonds = vec![];
    for (index2, order) in others {
        bonds.push((index1, index2, Bond::single()));
    }

    bonds
}

named!(gjf_connect_line<&str, Vec<(usize, usize, Bond)>>, do_parse!(
             opt!(space)           >>
    n:       unsigned_digit        >>
    others:  many0!(gjf_bond_pair) >>
             opt!(space)           >>
             line_ending           >>
    (
        build_bonds(n, others)
    )
));

#[test]
fn test_gjf_connectivity() {
    let (_, x) = gjf_connect_line(" 1,2 1.0 3 1.0 4 1.0 5 1.0\n").expect("gjf connectivity");
    assert_eq!(4, x.len());
}
// d03ec7e2-6cc0-475f-8fbc-d140db9ee4b2 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::2357368c-ab7e-4eb3-96a8-a8e4aba19bac][2357368c-ab7e-4eb3-96a8-a8e4aba19bac]]
named!(get_molecule_from<&str, Molecule>, do_parse!(
    link0: opt!(complete!(link0_section)) >>
    // route card
    route: route_section                  >>
    // molecule title
    title: title_section                  >>
    // charges and spin multipies
    chsps: charge_and_spin_line           >>
    // atom specification section
    atoms: many1!(gjf_atom_line)          >>
           blank_line                     >>
    // FIXME: how about gaussian extra input
    // connectivity section
    bonds: opt!(complete!(many0!(gjf_connect_line))) >>
    (
        {
            let mut mol = Molecule::new(title.trim());

            let mut map = HashMap::new();
            let mut i = 1;
            for a in atoms {
                let a = Atom::new(a.element_label, a.position);
                let n = mol.add_atom(a);
                map.insert(i, n);
                i += 1;
            }

            if let Some(bonds) = bonds {
                // FIXME: ugly
                for bonds in bonds {
                    for (i, j, b) in bonds {
                        // TODO: handle possible panic
                        let ai = map[&i];
                        let aj = map[&j];
                        mol.add_bond(ai, aj, b);
                    }
                }
            }

            mol
        }
    )
));

#[test]
fn test_gaussian_input() {
    let txt = "%chk=C5H12.chk
%nproc=8
%mem=5GB
#p opt freq=noraman nosymm B3LYP/6-31+G** test geom=connect

Title Card
Required

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

    let (_, mol) = get_molecule_from(txt).expect("gjf molecule");
    assert_eq!(17, mol.natoms());
    assert_eq!(16, mol.nbonds());
}

// TODO: atom properties
fn format_atom(a: &Atom) -> String {
    let [x, y, z] = a.position();

    format!(" {symbol:15} {x:14.8} {y:14.8} {z:14.8}\n",
            symbol=a.symbol(),
            x = x,
            y = y,
            z = z,
    )
}

// string representation in gaussian input file format
fn format_molecule(mol: &Molecule) -> String {
    let mut lines = String::new();

    let link0 = "%nproc=1\n%mem=20MW";
    let route = "#p sp scf=tight HF/3-21G* geom=connect test";
    lines.push_str(&format!("{}\n{}\n", link0, route));
    lines.push_str("\n");

    // title section
    lines.push_str("Title Card Required\n");
    lines.push_str("\n");

    // TODO: take from molecule
    lines.push_str("0 1\n");
    for a in mol.atoms() {
        let line = format_atom(&a);
        lines.push_str(&line);
    }

    // crystal vectors
    if let Some(lattice) = mol.lattice {
        let va = lattice.vector_a();
        let vb = lattice.vector_b();
        let vc = lattice.vector_c();
        for [x, y, z] in [va, vb, vc].iter() {
            lines.push_str(&format!(" TV              {:14.8}{:14.8}{:14.8}\n", x, y, z));
        }
    }

    // connectivity
    lines.push_str("\n");
    let mut map = HashMap::new();
    for (i, j, b) in mol.view_bonds() {
        let mut neighbors = map.entry(i).or_insert(vec![]);
        neighbors.push((j, b.order()));
    }
    for (i, a) in mol.view_atoms() {
        let mut line = format!("{:<5}", i);
        if let Some(neighbors) = map.get(&i) {
            for (j, o) in neighbors {
                line.push_str(&format!(" {:<5} {:<.1}", j, o));
            }
        }
        lines.push_str(&format!("{}\n", line));
    }

    lines.push_str("\n");
    lines
}
// 2357368c-ab7e-4eb3-96a8-a8e4aba19bac ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ac025fea-3d20-45f4-97eb-2969138a4716][ac025fea-3d20-45f4-97eb-2969138a4716]]
/// plain xyz coordinates with atom symbols
pub struct GaussInputFile();

/// References
/// http://gaussian.com/input/?tabid=0
impl ChemFileLike for GaussInputFile {
    fn ftype(&self) -> &str {
        "gaussian/input"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".gjf", ".com", ".gau"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule_from(chunk)
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        Ok(format_molecule(mol))
    }

    fn format(&self, mols: &Vec<Molecule>) -> Result<String> {
        let mut ms = vec![];
        for mol in mols {
            let m = format_molecule(mol);
            ms.push(m);
        }

        Ok(ms.join("--Link1--\n"))
    }
}
// ac025fea-3d20-45f4-97eb-2969138a4716 ends here
