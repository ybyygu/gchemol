// base

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*base][base:1]]
use super::*;
// base:1 ends here

// cell/elements

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*cell/elements][cell/elements:1]]
named!(poscar_cell_vectors<&str, [[f64; 3]; 3]>, sp!(do_parse!(
    va: xyz_array   >> eol >>
    vb: xyz_array   >> eol >>
    vc: xyz_array   >> eol >>
    (
        [va, vb, vc]
    )
)));

#[test]
fn test_poscar_cell_vectors() {
    let lines =
" 21.23300000  0.00000000  0.00000000
  0.00000000 26.60400000  0.00000000
  0.00000000  0.00000000 12.67600000
";

    let (_, x) = poscar_cell_vectors(lines).expect("POSCAR cell vectors");
    assert_eq!(21.233, x[0][0]);
}

named!(poscar_ion_types<&str, (Vec<&str>, Vec<usize>)>, sp!(do_parse!(
    syms: many1!(alpha)          >> eol >>
    nums: many1!(unsigned_digit) >> eol >>
    ((syms, nums))
)));

#[test]
fn test_formats_vasp_poscar_ion_types() {
    let lines = " O    Si   C    N    H
 225  112   8    1    19 \n";
    let (_, v) = poscar_ion_types(lines).expect("POSCAR ion types");
    assert_eq!(5, v.0.len());
}
// cell/elements:1 ends here

// coordinates

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*coordinates][coordinates:1]]
// Selective dynamics -- optional, can be omitted
// only the first character is relevant
named!(selective_dynamics<&str, &str>, do_parse!(
    s: tag_no_case!("S") >> read_line >>
    (s)
));

// Direct/Cartesian -- lattice coordinates type
// only the first character is relevant
named!(direct_or_catersian<&str, &str>, do_parse!(
    d: alt!(tag_no_case!("D") | tag_no_case!("C")) >> read_line >>
    (d)
));

// combine two parsers
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

    let (_, (s, d)) = poscar_select_direct(lines).expect("poscar selective/direct");
    assert_eq!(true, s);
    assert_eq!(true, d);

    let (_, (s, d)) = poscar_select_direct("Direct\n").expect("poscar direct");
    assert_eq!(false, s);
    assert_eq!(true, d);
    let (_, (s, d)) = poscar_select_direct("Cartesian\n").expect("poscar catersian");
    assert_eq!(false, s);
    assert_eq!(false, d);
}

/// Consume three chars in selective dynamics flag (T/F) separated by one or more spaces
/// Return the frozen flag array
named!(pub selective_dynamics_flags<&str, [bool; 3]>, do_parse!(
    x: one_of!("TF") >> space >>
    y: one_of!("TF") >> space >>
    z: one_of!("TF") >>
    (
        [x == 'T', y == 'T', z == 'T']
    )
));

#[test]
fn test_poscar_select_dynamics() {
    let (_, x) = selective_dynamics_flags("T T F").unwrap();
    assert_eq!(x, [true, true, false]);
}

// 0.05185     0.39121     0.29921  T T T # O
// 0.81339     0.57337     0.68777  T T T # O
named!(poscar_position<&str, ([f64; 3], Option<[bool; 3]>)>, sp!(do_parse!(
    position: xyz_array                                   >>
    //frozen  : opt!(complete!(selective_dynamics_flags)) >> read_line >>
    frozen  : opt!(selective_dynamics_flags) >> read_line >>
    ((position, frozen))
)));

#[test]
fn test_poscar_position() {
    let line = "     0.05185     0.39121     0.29921  T T T # O \n";
    let (_, (position, sflags)) = poscar_position(line).expect("POSCAR position style 1");
    assert_eq!(0.05185, position[0]);
    assert_eq!(0.39121, position[1]);
    assert_eq!(0.29921, position[2]);
    assert_eq!(Some([true, true, true]), sflags);

    let line = "     0.05185     0.39121     0.29921\n";
    let (_, (position, sflags)) = poscar_position(line).expect("POSCAR position style 1");
    assert_eq!(None, sflags);
}
// coordinates:1 ends here

// parse molecule

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*parse%20molecule][parse molecule:1]]
/// Read Molecule from stream in VASP/POSCAR format
named!(get_molecule_from<&str, Molecule>, do_parse!(
    title            : sp!(read_line)          >>
    lattice_constant : read_f64                >>
    cell_vectors     : poscar_cell_vectors     >>
    ion_types        : poscar_ion_types        >>
    select_direct    : poscar_select_direct    >>
    ion_positions    : many1!(poscar_position) >>
    (
        {
            let selective_dynamics = select_direct.0;
            let direct_coordinates = select_direct.1;

            let mut mol = Molecule::new(title);
            let mut lat = Lattice::new(cell_vectors);
            lat.scale_by(lattice_constant);

            let mut symbols = vec![];
            let (syms, nums) = ion_types;

            for (&sym, &num) in syms.iter().zip(nums.iter()) {
                for _ in 0..num {
                    symbols.push(sym);
                }
            }

            if symbols.len() != ion_positions.len() {
                eprintln!("WARNING: some ions data not correctly parsed!");
            }

            for (&sym, (pos, sflags)) in symbols.iter().zip(ion_positions) {
                let mut pos = pos;
                if direct_coordinates {
                    pos = lat.to_cart(pos);
                }
                let mut a = Atom::new(sym, pos);
                // FIXME: just a temporary workaround
                if sflags.is_some() {
                    a.properties.store(POSCAR_SFLAGS_KEY, sflags.unwrap());
                };

                mol.add_atom(a);
            }
            mol.set_lattice(lat);

            mol
        }
    )
));

#[test]
fn test_poscar_molecule() {
    let lines =
"title
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
// parse molecule:1 ends here

// format molecule

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*format%20molecule][format molecule:1]]
const POSCAR_SFLAGS_KEY: &str = "vasp/poscar/sflags";

fn format_molecule(mol: &Molecule) -> String {
    let mut lines = String::new();
    let title = if mol.name.is_empty() {"generated by gchemol"} else {&mol.name};
    lines.push_str(&format!("{}\n", title));
    lines.push_str("1.0\n");
    let mut lattice = mol.lattice.expect("poscar lattice");
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

    // write fractional coordinates for improving accuracy
    lines.push_str("Selective dynamics\nDirect\n");
    for a in mol.atoms() {
        let [x, y, z] = lattice.to_frac(a.position());
        // FIXME: just a temporary workaround
        let line = if a.properties.contains_key(POSCAR_SFLAGS_KEY) {
            let sflags: [bool; 3] = a.properties.load(POSCAR_SFLAGS_KEY).expect("vasp selective_dynamics flags");
            format!("{x:18.12} {y:18.12} {z:18.12} {fx} {fy} {fz}\n",
                    x=x,
                    y=y,
                    z=z,
                    fx= if sflags[0] {"T"} else {"F"},
                    fy= if sflags[1] {"T"} else {"F"},
                    fz= if sflags[2] {"T"} else {"F"},
            )
        } else {
            format!("{x:18.12} {y:18.12} {z:18.12} T T T\n",
                    x=x,
                    y=y,
                    z=z)
        };

        lines.push_str(&line);
    }

    // final blank line
    lines.push_str("\n");
    // TODO: write velocities

    lines
}

// Panic if symbols is empty
fn count_symbols(symbols: Vec<&str>) -> Vec<(&str, usize)> {
    let mut lines = String::new();

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
// format molecule:1 ends here

// chemfile

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*chemfile][chemfile:1]]
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
    let mols = poscar.parse(Path::new("tests/files/vasp/POSCAR")).unwrap();
    assert_eq!(1, mols.len());
    assert_eq!(365, mols[0].natoms());
}
// chemfile:1 ends here
