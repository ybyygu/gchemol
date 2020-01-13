// header

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*header][header:1]]
// Data will be parsed:
// Lattice, Atoms, Bonds
// header:1 ends here

// base

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*base][base:1]]
use super::*;
// base:1 ends here

// cell loop

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*cell loop][cell loop:1]]
use nom::recognize_float;

/// Recognizes a float point number with uncertainty brackets
/// # Example
/// 10.154(2)
named!(double_cif<&str, f64>, do_parse!(
    v: recognize_float                                 >>
    o: opt!(delimited!(char!('('), digit, char!(')'))) >>

    (
        {
            let s = if let Some(o) = o {
                v.to_owned() + o
            } else {
                v.to_owned()
            };

            s.parse().expect("cif uncertainty number")
        }
    )
));

#[test]
fn test_cif_float_number() {
    let (_, v) = double_cif("0.3916\n").expect("cif float1");
    assert_eq!(v, 0.3916);
    let (_, v) = double_cif("0.391(6)\n").expect("cif float2");
    assert_eq!(v, 0.3916);
}

/// Recognizes a float value with a preceeding tag
fn tagged_f64<'a>(input: &'a str, tag: &'a str) -> IResult<&'a str, f64> {
    sp!(input, preceded!(tag!(tag), double_cif))
}

#[test]
fn test_tagged_f64() {
    let (_, v) = tagged_f64(" abc 4.1 \n", "abc").expect("cif tagged f64");
    assert_eq!(4.1, v);
}

named!(cell_params<&str, (f64, f64, f64, f64, f64, f64)>, ws!(permutation!(
    call!(tagged_f64, "_cell_length_a"),
    call!(tagged_f64, "_cell_length_b"),
    call!(tagged_f64, "_cell_length_c"),
    call!(tagged_f64, "_cell_angle_alpha"),
    call!(tagged_f64, "_cell_angle_beta"),
    call!(tagged_f64, "_cell_angle_gamma")
)));

/// Read crystal cell
named!(read_cell<&str, Lattice>, do_parse!(
    params: cell_params >>
    (
        Lattice::from_params(
            params.0,
            params.1,
            params.2,
            params.3,
            params.4,
            params.5,
        )
    )
));

#[test]
fn test_cif_cell_loop() {
    // Allow data in random order and with blank line
    let lines = "_cell_length_a                    18.094(0)
_cell_length_c                    7.5240

_cell_length_b                    20.5160
_cell_angle_alpha                 90.0000
_cell_angle_beta                  90.0000
_cell_angle_gamma                 90.0000
loop_
";

    let (_, param) = cell_params(lines).expect("cif cell");
    assert_eq!(param.1, 20.5160);
}
// cell loop:1 ends here

// atom sites loop
// # Example
// loop_
// _atom_site_label
// _atom_site_type_symbol
// _atom_site_fract_x
// _atom_site_fract_y
// _atom_site_fract_z
// Cu1 Cu 0.20761(4) 0.65105(3) 0.41306(4)
// O1 O 0.4125(2) 0.6749(2) 0.5651(3)
// O2 O 0.1662(2) 0.4540(2) 0.3821(3)
// O3 O 0.4141(4) 0.3916(3) 0.6360(4)
// N1 N 0.2759(3) 0.8588(2) 0.4883(3)
// ...

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*atom sites loop][atom sites loop:1]]
named!(atom_site_header<&str, &str>, preceded!(
    tag!("_atom_site_"),
    not_space
));

// ["label", "type_symbol", "fract_x", "fract_y", "fract_z", "U_iso_or_equiv", "adp_type", "occupancy"]
named!(atom_site_headers<&str, Vec<&str>>, do_parse!(
    headers: many1!(ws!(atom_site_header)) >>
    (
        headers
    )
));

#[test]
fn test_cif_site_headers() {
    let lines = "\
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y

_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy
Si1    Si    0.30070   0.07240   0.04120   0.00000  Uiso   1.00 \n";

    let (_, h) = atom_site_headers(lines).expect("cif atom site headers");
    assert_eq!(h.len(), 8);
}

/// Read atoms in cif _atom_site loop
fn read_atoms<'a>(input: &'a str) -> IResult<&'a str, Vec<Atom>> {
    let (rest, headers) = atom_site_headers(input)?;
    if headers.len() <= 4 {
        println!("{:?}", rest);
        eprintln!("cif formats: not enough columns in atom site loop");
        // return Err(nom_failure!(rest));
        panic!("nom failure")
    }

    // column header loopup table
    // Example
    // -------
    //   0        1         2       3      4            5          6         7
    // label type_symbol fract_x fract_y fract_z U_iso_or_equiv adp_type occupancy
    let table: HashMap<_, _> = HashMap::from(headers.iter().zip(0..).collect());
    let ifx = *table.get(&"fract_x").expect("fract x col");
    let ify = *table.get(&"fract_y").expect("fract y col");
    let ifz = *table.get(&"fract_z").expect("fract z col");
    // column index to atom label
    let ilbl = *table.get(&"label").expect("atom label col");
    // TODO: column index to element symbol, which is optional
    let isym = *table.get(&"type_symbol").expect("atom symbol col");

    do_parse!(
        rest,
        rows: many1!(terminated!(
            count!(sp!(not_space), headers.len()),
            sp!(line_ending)
        )) >> ({
            let mut atoms = vec![];
            for row in rows {
                let line = format!("{} {} {}\n", row[ifx], row[ify], row[ifz]);
                let (_, (fx, fy, fz)) = sp!(&line, tuple!(double_cif, double_cif, double_cif))
                    .map_err(|e| {
                        eprintln!("failed to parse coordinates: {:?}", e);
                        panic!("no failure");
                    })?;

                let lbl = row[ilbl];
                let sym = row[isym];
                // TODO: assign atom label
                let a = Atom::build().symbol(sym).position(fx, fy, fz).finish();
                atoms.push(a);
            }

            atoms
        })
    )
}

#[test]
fn test_cif_atoms() {
    let lines = "_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_U_iso_or_equiv
_atom_site_adp_type
_atom_site_occupancy
Si1    Si    0.30070   0.07240   0.04120   0.00000  Uiso   1.00
Si2    Si    0.30370   0.30880   0.04610   0.00000  Uiso   1.00
O3     O     0.12430   0.41700   0.42870   0.00000  Uiso   1.00
O4     O     0.12260   0.19540   0.42540   0.00000  Uiso   1.00
O5     O     0.23620   0.12240   0.98650   0.00000  Uiso   1.00
Si6    Si    0.80070   0.57240   0.04120   0.00000  Uiso   1.00
Si7    Si    0.80370   0.80880   0.04610   0.00000  Uiso   1.00
O8     O     0.62430   0.91700   0.42870   0.00000  Uiso   1.00
O9     O     0.62260   0.69540   0.42540   0.00000  Uiso   1.00
O10    O     0.73620   0.62240   0.98650   0.00000  Uiso   1.00
Si11   Si    0.69930   0.92760   0.54120   0.00000  Uiso   1.00
Si12   Si    0.69630   0.69120   0.54610   0.00000  Uiso   1.00
\n";
    let (r, v) = read_atoms(lines).expect("cif atom site loop");
    assert_eq!(12, v.len());
}
// atom sites loop:1 ends here

// bond loop
// # Example
// loop_
// _geom_bond_atom_site_label_1
// _geom_bond_atom_site_label_2
// _geom_bond_distance
// _geom_bond_site_symmetry_2
// _ccdc_geom_bond_type
// Si1    O140    1.629   .     S
// Si1    O128    1.624   .     S
// Si1    O5      1.607   1_554 S
// Si1    O18     1.614   1_554 S
// Si2    O86     1.587   .     S
// ...

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*bond loop][bond loop:1]]
named!(geom_bond_header<&str, &str>, preceded!(
    alt!(
        tag!("_geom_bond_") |
        tag!("_ccdc_geom_bond_")
    ),
    not_space
));

named!(geom_bond_headers<&str, Vec<&str>>, ws!(preceded!(
    tag!("loop_"),
    many1!(geom_bond_header)
)));

#[test]
fn test_cif_bond_header() {
    let txt = "loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_geom_bond_distance
_geom_bond_site_symmetry_2
_ccdc_geom_bond_type
# END
";

    let (_, x) = geom_bond_headers(txt).expect("cif bond headers");
    assert_eq!(5, x.len());
}
// bond loop:1 ends here

// molecule

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*molecule][molecule:1]]
// The first line
named!(cif_title<&str, &str>, preceded!(
    tag!("data_"),
    not_space
));

named!(cell_start<&str, &str>, sp!(
    tag!("_cell_length_")
));

named!(atom_loop_start<&str, &str>, sp!(
    tag!("_atom_site_")
));


/// Create Molecule object from cif stream
fn read_molecule(input: &str) -> IResult<&str, Molecule> {
    // goto block start
    let (input, _) = read_many_until(input, cif_title)?;
    let (input, title) = cif_title(input)?;

    // read cell parameters
    let (input, _) = read_many_until(input, cell_start)?;
    let (input, lat) = read_cell(input)?;

    // read atoms
    let (input, _) = read_many_until(input, atom_loop_start)?;
    let (input, atoms) = read_atoms(input)?;

    // create molecule
    let mut mol = Molecule::new(title);
    for mut a in atoms {
        let p = lat.to_cart(a.position());
        a.set_position(p);
        mol.add_atom(a);
    }

    mol.set_lattice(lat);

    // TODO: add bonds

    Ok((input, mol))
}
// molecule:1 ends here

// chemfile

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*chemfile][chemfile:1]]
pub struct CifFile();

impl ChemFileLike for CifFile {
    fn ftype(&self) -> &str {
        "text/cif"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".cif"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        read_molecule(chunk)
    }

    /// Represent molecule in .cif format
    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        let mut lines = String::new();

        // 1. meta inforation
        lines.push_str("data_test\n");
        lines.push_str("_audit_creation_method            'gchemol'\n");
        lines.push_str("_symmetry_space_group_name_H-M    'P1'\n");
        lines.push_str("_symmetry_Int_Tables_number       1\n");
        lines.push_str("_symmetry_cell_setting            triclinic\n");
        lines.push_str("\n");

        // 2. cell parameters
        lines.push_str("loop_\n");
        lines.push_str("_symmetry_equiv_pos_as_xyz\n");
        lines.push_str(" x,y,z\n");

        let mut lat = mol.lattice.ok_or(format_err!("Not a periodic moelcule."))?;
        let [a, b, c] = lat.lengths();
        let [alpha, beta, gamma] = lat.angles();
        lines.push_str(&format!("_cell_length_a     {:10.4}\n", a));
        lines.push_str(&format!("_cell_length_b     {:10.4}\n", b));
        lines.push_str(&format!("_cell_length_c     {:10.4}\n", c));
        lines.push_str(&format!("_cell_angle_alpha  {:10.4}\n", alpha));
        lines.push_str(&format!("_cell_angle_beta   {:10.4}\n", beta));
        lines.push_str(&format!("_cell_angle_gamma  {:10.4}\n", gamma));
        lines.push_str("\n");

        // 3. atom fractional coordinates
        lines.push_str("loop_\n");
        lines.push_str("_atom_site_type_symbol\n");
        lines.push_str("_atom_site_label\n");
        lines.push_str("_atom_site_fract_x\n");
        lines.push_str("_atom_site_fract_y\n");
        lines.push_str("_atom_site_fract_z\n");

        for a in mol.atoms() {
            let position = a.position();
            let symbol = a.symbol();
            let name = a.label();
            let [fx, fy, fz] = lat.to_frac(position);
            let s = format!("{:4}{:6}{:12.5}{:12.5}{:12.5}\n",
                            symbol,
                            name,
                            fx,
                            fy,
                            fz);
            lines.push_str(&s);
        }

        // 4. bonds
        // if mol.nbonds() > 0 {
        //     lines.push_str("loop_\n");
        //     lines.push_str("_geom_bond_atom_site_label_1\n");
        //     lines.push_str("_geom_bond_atom_site_label_2\n");
        //     lines.push_str("_geom_bond_distance\n");
        //     lines.push_str("_geom_bond_site_symmetry_2\n");
        //     lines.push_str("_ccdc_geom_bond_type\n");
        //     for bond in mol.bonds() {
        //         let symbol1 = frame.symbols.get(&current).unwrap();
        //         let name1 = format!("{}{}", symbol1, current);
        //         let p1 = frame.positions.get(&current).unwrap();
        //         let p1 = Point3::new(p1[0], p1[1], p1[2]) - cell_origin;

        //         let connected = frame.neighbors.get(&current).unwrap();
        //         for other in connected {
        //             if *other > current {
        //                 let symbol2 = frame.symbols.get(&other).unwrap();
        //                 let name2 = format!("{}{}", symbol2, other);
        //                 let p2 = frame.positions.get(&other).unwrap();
        //                 let p2 = Point3::new(p2[0], p2[1], p2[2]) - cell_origin;
        //                 let (image, distance) = get_nearest_image(cell, p1, p2);
        //                 if image.x == 0. && image.y == 0. && image.z == 0. {
        //                     lines.push_str(&format!("{:6} {:6} {:6.3} {:6} S\n", name1, name2, distance, "."));
        //                 } else {
        //                     let symcode = get_image_symcode(image);
        //                     lines.push_str(&format!("{:6} {:6} {:6.3} {:6} S\n", name1, name2, distance, symcode));
        //                     let (image, distance) = get_nearest_image(cell, p2, p1);
        //                     let symcode = get_image_symcode(image);
        //                     lines.push_str(&format!("{:6} {:6} {:6.3} {:6} S\n", name2, name1, distance, symcode));
        //                 }
        //             }
        //         }
        //     }
        // }

        Ok(lines)
    }
}
// chemfile:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*test][test:1]]
#[test]
fn test_formats_cif() {
    let file = CifFile();
    let fname = "tests/files/cif/babel.cif";
    let mols = file.parse(Path::new(fname)).expect("babel.cif");
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].natoms(), 34);

    let fname = "tests/files/cif/MS-MOR.cif";
    let mols = file.parse(Path::new(fname)).expect("MS-MOR.cif");
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].natoms(), 144);

    let fname = "tests/files/cif/IZA-LTL.cif";
    let mols = file.parse(Path::new(fname)).expect("cif IZA");
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].natoms(), 8);

    let fname = "tests/files/cif/ccdc.cif";
    let mols = file.parse(Path::new(fname)).expect("cif ccdc");
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].natoms(), 41);

    let fname = "tests/files/cif/quinone.cif";
    let mols = file.parse(Path::new(fname)).expect("cif quinone");
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].natoms(), 16);
}
// test:1 ends here
