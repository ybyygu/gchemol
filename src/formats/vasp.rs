// [[file:~/Workspace/Programming/gchemol/formats.note::580d4de7-5923-4885-87d5-20b24d8a703b][580d4de7-5923-4885-87d5-20b24d8a703b]]
// VASP POSCAR file format:
// http://cms.mpi.univie.ac.at/vasp/guide/node59.html
// 580d4de7-5923-4885-87d5-20b24d8a703b ends here

// [[file:~/Workspace/Programming/gchemol/formats.note::70e8dd41-45af-4712-b095-802079ac6eb4][70e8dd41-45af-4712-b095-802079ac6eb4]]
use super::*;

named!(poscar_lattice_constant<&str, f64>, terminated!(
    sp!(double_s),
    sp!(line_ending)
));

#[test]
fn test_poscar_lattice_constant() {
    let (_, x) = poscar_lattice_constant(" 1.0000 \n")
        .expect("POSCAR lattice constant");
    assert_eq!(1.000, x);
}

named!(poscar_cell_vectors<&str, [[f64; 3]; 3]>, do_parse!(
    va: sp!(xyz_array)   >>
        sp!(line_ending) >>
    vb: sp!(xyz_array)   >>
        sp!(line_ending) >>
    vc: sp!(xyz_array)   >>
        sp!(line_ending) >>
    (
        [va, vb, vc]
    )
));

#[test]
fn test_poscar_cell_vectors() {
    let lines = " 21.23300000  0.00000000  0.00000000
  0.00000000 26.60400000  0.00000000
  0.00000000  0.00000000 12.67600000 \n";
    let (_, x) = poscar_cell_vectors(lines)
        .expect("POSCAR cell vectors");
    assert_eq!(21.233, x[0][0]);
}

named!(poscar_ion_types<&str, (Vec<&str>, Vec<usize>)>, do_parse!(
    syms: many1!(sp!(alpha))          >>
    sp!(line_ending)                       >>
    nums: many1!(sp!(unsigned_digit)) >>
    sp!(line_ending)                       >>
    ((syms, nums))
));

#[test]
fn test_formats_vasp_poscar_ion_types() {
    let lines = " O    Si   C    N    H
 225  112   8    1    19 \n";
    let (r, v) = poscar_ion_types(lines)
        .expect("POSCAR ion types");
    assert_eq!(5, v.0.len());
}

// Selective dynamics -- optional, can be omitted
// only the first character is relevant
named!(selective_dynamics<&str, &str>, do_parse!(
    s: tag_no_case!("S") >> read_until_eol >>
    (s)
));

// Direct/Cartesian -- lattice coordinates type
// only the first character is relevant
named!(direct_or_catersian<&str, &str>, do_parse!(
    d: sp!(alt!(tag_no_case!("D") | tag_no_case!("C"))) >> read_until_eol >>
    (d)
));

// combine two above parsers
named!(poscar_select_direct<&str, (bool, bool)>, do_parse!(
    s: opt!(complete!(selective_dynamics)) >>
    d: direct_or_catersian >>
    (
        {
            let d = d.to_lowercase();
            (s.is_some(), d == "d")
        }
    )
));

#[test]
fn test_poscar_select_direct() {
    let lines = "Selective dynamics
Direct\n";

    let (_, (s, d)) = poscar_select_direct(lines)
        .expect("poscar selective/direct");
    assert_eq!(true, s);
    assert_eq!(true, d);

    let (_, (s, d)) = poscar_select_direct("Direct\n").expect("poscar direct");
    assert_eq!(false, s);
    assert_eq!(true, d);
    let (_, (s, d)) = poscar_select_direct("Cartesian\n").expect("poscar catersian");
    assert_eq!(false, s);
    assert_eq!(false, d);
}

// TODO: selective dynamics labels
// 0.05185     0.39121     0.29921  T T T # O
// 0.81339     0.57337     0.68777  T T T # O
named!(poscar_position<&str, ([f64; 3], &str)>, do_parse!(
    position: sp!(xyz_array)                    >>
    remained: read_until_eol                    >>
    ((position, remained))
));

#[test]
fn test_poscar_position() {
    let x = poscar_position("     0.05185     0.39121     0.29921  T T T # O \n")
        .expect("POSCAR position style 1");
    let x = poscar_position("     0.05185     0.39121     0.29921 \n")
        .expect("POSCAR position style 1");
}

// TODO: read velocities
named!(get_molecule_from<&str, Molecule>, do_parse!(
    title            : sp!(read_until_eol)        >>
    lattice_constant : poscar_lattice_constant    >>
    cell_vectors     : poscar_cell_vectors        >>
    ion_types        : poscar_ion_types           >>
    select_direct    : poscar_select_direct       >>
    ion_positions    : many1!(poscar_position)    >>
    (
        {
            let selective_dynamics = select_direct.0;
            let direct_coordinates = select_direct.1;
            let mut mol = Molecule::new(title);
            let mut tvs = cell_vectors;
            let mut lat = Lattice::new(tvs);
            lat.scale_by(lattice_constant);

            let mut symbols = vec![];
            let (syms, nums) = ion_types;

            for (&sym, &num) in syms.iter().zip(nums.iter()) {
                for _ in 0..num {
                    symbols.push(sym);
                }
            }

            if symbols.len() != ion_positions.len() {
                eprintln!("Inconsistency: some ions data not correctly parsed.");
            }

            for (&sym, (pos, codes)) in symbols.iter().zip(ion_positions) {
                let mut pos = pos;
                if direct_coordinates {
                    pos = lat.to_cart(pos);
                }
                let a = Atom::new(sym, pos);
                mol.add_atom(a);
            }
            mol.set_lattice(lat);

            mol
        }
    )
));

#[test]
fn test_poscar_molecule() {
    let lines = "title
1.0
 21.23300000  0.00000000  0.00000000
  0.00000000 26.60400000  0.00000000
  0.00000000  0.00000000 12.67600000
 O    Si   C    N    H
225  112   8    1    19
Selective dynamics
Direct
     0.05185     0.39121     0.29921  T T T # O
     0.81339     0.57337     0.68777  T T T # O
     0.73422     0.23229     0.85313  T T T # O
     0.02246     0.05156     0.49349  T T T # O
     0.64451     0.66726     0.17130  T T T # O
     0.05185     0.07337     0.29921  T T T # O
     0.60095     0.57471     0.17096  T T T # O
     0.64451     0.66726     0.81569  T T T # O
     0.33416     0.64745     0.88951  T T T # O
     0.33416     0.31713     0.09747  T T T # O
     0.93262     0.92263     0.99349  T T T # O
     0.43262     0.79195     0.99349  T T T # O
     0.73422     0.73229     0.13386  T T T # O
     0.22073     0.66726     0.81569  T T T # O
\n";

    let (_, mol) = get_molecule_from(lines).expect("poscar molecule");
    assert_eq!(14, mol.natoms());
    assert!(mol.lattice.is_some());
}

// Panic if symbols is empty
fn count_symbols(symbols: Vec<&str>) -> Vec<(&str, usize)> {
    use indexmap::IndexMap;

    let mut lines = String::new();

    // let mut counting_symbols = IndexMap::new();
    // for sym in mol.symbols() {
    //     *counting_symbols.entry(sym).or_insert(0) += 1;
    // }

    let mut syms1 = symbols.iter();
    let mut syms2 = symbols.iter().skip(1);
    let mut counts = vec![];

    let mut c = 1;
    let mut s = symbols[0];
    for (&sym1, &sym2) in syms1.zip(syms2) {
        if sym2 == sym1 {
            c += 1;
        } else {
            counts.push((sym1, c));
            c = 1;
        }
        s = sym2;
    }
    // append the last piece
    counts.push((s, c));

    counts
}

#[test]
fn test_poscar_symbols_counts() {
    let symbols = ["C", "C", "C", "H", "O", "O", "C"];
    let x = count_symbols(symbols.to_vec());
    assert_eq!([("C", 3), ("H", 1), ("O", 2), ("C", 1)].to_vec(), x);

    let symbols = ["C", "C"];
    let x = count_symbols(symbols.to_vec());
    assert_eq!([("C", 2)].to_vec(), x);

    let symbols = ["C", "H"];
    let x = count_symbols(symbols.to_vec());
    assert_eq!([("C", 1), ("H", 1)].to_vec(), x);

    let symbols = ["C"];
    let x = count_symbols(symbols.to_vec());
    assert_eq!([("C", 1)].to_vec(), x);
}

fn format_molecule(mol: &Molecule) -> String {
    let mut lines = String::new();
    let title = if mol.name.is_empty() {"generated by gchemol"} else {&mol.name};
    lines.push_str(&format!("{}\n", title));
    lines.push_str("1.0\n");
    let lattice = mol.lattice.expect("poscar lattice");
    let va = lattice.vector_a();
    let vb = lattice.vector_b();
    let vc = lattice.vector_c();

    for v in [va, vb, vc].into_iter() {
        let line = format!("{:12.8}{:12.8}{:12.8}\n", v[0], v[1], v[2]);
        lines.push_str(&line);
    }

    // atom symbols and counts
    let mut line1 = String::new();
    let mut line2 = String::new();
    for (s, n) in count_symbols(mol.symbols()) {
        line1.push_str(&format!(" {:^4}", s));
        line2.push_str(&format!(" {:^4}", n));
    }
    lines.push_str(&format!("{}\n", line1));
    lines.push_str(&format!("{}\n", line2));

    // TODO: write fractional coordinates?
    lines.push_str("Selective dynamics\nCatersian\n");
    for a in mol.atoms() {
        let [x, y, z] = a.position();
        // TODO
        let line = format!("{x:12.5}{y:12.5}{z:12.5} T T T\n",
                           x=x,
                           y=y,
                           z=z);
        lines.push_str(&line);
    }

    // final blank line
    lines.push_str("\n");
    // TODO: write velocities

    lines
}
// 70e8dd41-45af-4712-b095-802079ac6eb4 ends here

// [[file:~/Workspace/Programming/gchemol/formats.note::e71c5ee4-27e3-4c03-b626-e7ae375b4510][e71c5ee4-27e3-4c03-b626-e7ae375b4510]]
use io;

pub struct PoscarFile();

impl ChemFileLike for PoscarFile {
    fn ftype(&self) -> &str {
        "vasp/poscar"
    }

    fn extensions(&self) -> Vec<&str> {
        vec!["POSCAR", "CONTCAR", ".poscar", ".vasp"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule_from(chunk)
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        Ok(format_molecule(mol))
    }
}

#[test]
fn test_vasp_poscar_parse() {
    let poscar = PoscarFile();
    let mols = poscar.parse("tests/files/vasp/POSCAR").unwrap();
    assert_eq!(1, mols.len());
    assert_eq!(365, mols[0].natoms());
}
// e71c5ee4-27e3-4c03-b626-e7ae375b4510 ends here
