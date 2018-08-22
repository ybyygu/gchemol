// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::580d4de7-5923-4885-87d5-20b24d8a703b][580d4de7-5923-4885-87d5-20b24d8a703b]]
// DFTBPlus GEN format
//
// a generalized xyz format for describing bare clusters and periodic structures
//
// For details, see appendix C of offical user manual
// http://www.dftbplus.org/fileadmin/DFTBPLUS/public/dftbplus/latest/manual.pdf
// 580d4de7-5923-4885-87d5-20b24d8a703b ends here

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::00b6b7ea-25e9-4e42-ba0a-2af4a8b63013][00b6b7ea-25e9-4e42-ba0a-2af4a8b63013]]
use nom;
use super::*;

// The first line of the file contains the number of atoms, n, followed by the
// type of geometry.
//
// C for cluster (non-periodic), S for supercell in Cartesian coordinates or F
// for supercell in fractions of the lattice vectors.
named!(get_geom_meta<&str, (usize, char)>, do_parse!(
    natoms    : sp!(unsigned_digit) >>
    geom_type : sp!(one_of!("CSF")) >>
                sp!(line_ending)    >>
    (
        natoms,
        geom_type
    )
));

#[test]
fn test_dftb_geom_meta() {
    let line = "12 C\n";
    let (_, (n, x)) = get_geom_meta(line).expect("dftb geom meta");
    assert_eq!(12, n);
}

// The second line contains the chemical symbols of the elements present
// separated by one or more spaces.
named!(get_element_symbols<&str, Vec<&str>>, do_parse!(
    symbols: separated_nonempty_list!(space, sp!(alpha)) >>
             sp!(line_ending) >>
    (symbols)
));

#[test]
fn test_dftb_element_symbols() {
    let line = " C  H O\n";
    let (_, x) = get_element_symbols(line).expect("dftb element symbols");
    assert_eq!(3, x.len());
}

// The following n lines contain a list of the atoms.
named!(get_atom_type_and_position<&str, (usize, [f64; 3])>, do_parse!(
    index: sp!(unsigned_digit) >>
    kind: sp!(unsigned_digit)  >>
    position: sp!(xyz_array)   >>
    sp!(line_ending)           >>
    (
        kind,
        position
    )
));

#[test]
fn test_dftb_atom_type_and_positions() {
    let line = "    1 1    0.1875714333E+02    0.1561236879E+02    0.7500000000E+01 \n";
    let (_, (k, p)) = get_atom_type_and_position(line).expect("dftb atom type/positions");
    assert_eq!(1, k);
    assert_relative_eq!(18.75714, p[0], epsilon=1e-3);
}

// crystal lattice data: cell origin and vectors
named!(get_lattice<&str, Lattice>, do_parse!(
    origin : sp!(xyz_array)   >>
             sp!(line_ending) >>
    tvx    : sp!(xyz_array)   >>
             sp!(line_ending) >>
    tvy    : sp!(xyz_array)   >>
             sp!(line_ending) >>
    tvz    : sp!(xyz_array)   >>
             sp!(line_ending) >>
    (
        {
            let mut lat = Lattice::new(
                [
                    tvx,
                    tvy,
                    tvz
                ]
            );
            lat.set_origin(origin);
            lat
        }
    )
));

#[test]
fn test_dftb_lattice() {
    let lines = "0.000000   0.000000   0.000000
2.713546   2.713546   0.000000
0.000000   2.713546   2.713546
2.713546   0.000000   2.713546\n";
    let (_, x) = get_lattice(lines).expect("dftb lattice");
    assert_relative_eq!(2.713546, x.vector_a()[0], epsilon=1e-3);
}
// 00b6b7ea-25e9-4e42-ba0a-2af4a8b63013 ends here

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::cb010521-689c-4b07-9cba-4e2dbea7fac3][cb010521-689c-4b07-9cba-4e2dbea7fac3]]
use indexmap::IndexSet;
use itertools::Itertools;

fn get_molecule(input: &str) -> IResult<&str, Molecule> {
    let (rest, (natoms, geom_type)) = get_geom_meta(input)?;

    // check geometry type
    let (periodic, fractional) = match geom_type {
        'C' => {
            debug!("dftb geom type: cluster");
            (false, false)
        },
        'S' => {
            debug!("dftb geom type: periodic + cartesian coordinates");
            (true, false)
        },
        'F' => {
            debug!("dftb geom type: periodic + fractional coordinates");
            (true, true)
        },
        _   => {
            error!("invalid geometry type");
            // FIXME: nom error
            return Err(nom::Err::Failure(error_position!(rest, nom::ErrorKind::Custom(29))));
        }
    };
    // read element types
    let (mut rest, symbols) = get_element_symbols(rest)?;
    let mut atoms = vec![];
    for _ in 0..natoms {
        let (r, (kind, position)) = get_atom_type_and_position(rest)?;
        if kind > symbols.len() {
            error!("element type index out of range");
            // FIXME: nom error
            return Err(nom::Err::Failure(error_position!(rest, nom::ErrorKind::Custom(30))));
        }
        let sym = symbols[kind - 1];
        let a = Atom::new(sym, position);
        atoms.push(a);
        rest = r;
    }

    let mut mol = Molecule::new("from dftb gen");
    // read lattice
    if periodic {
        let (r, lattice) = get_lattice(rest)?;
        // convert fractional coordinates to cartesian
        if fractional {
            for i in 0..natoms {
                let pos = atoms[i].position();
                let pos = lattice.to_cart(pos);
                atoms[i].set_position(pos);
            }
        }

        // update the remaing content
        rest = r;
        // assign lattice
        mol.set_lattice(lattice);
    }

    for a in atoms {
        mol.add_atom(a);
    }

    Ok((rest, mol))
}

#[test]
fn test_dftb_get_molecule() {
    let lines = read_file("tests/files/dftb/geo_end.gen").expect("dftb gen");
    let (_, x) = get_molecule(&lines).expect("dftb molecule");
    assert_eq!(12, x.natoms());
}

// format molecule in dftb gen-format
fn format_molecule(mol: &Molecule) -> Result<String> {
    let mut lines = String::new();
    // the first line
    let geom_type = if mol.lattice.is_none() {
        "C"
    } else {
        "S"
    };
    let line = format!("{} {}\n", mol.natoms(), geom_type);
    lines.push_str(&line);

    // the second line
    let symbols: IndexSet<_> = mol.symbols().into_iter().collect();
    let line = format!(
        "{}\n",
        symbols.iter().join(" ")
    );
    lines.push_str(&line);

    // the following n lines for atom positions
    for (i, a) in mol.view_atoms() {
        let position = a.position();
        let sym = a.symbol();
        let (tag, _) = symbols
            .get_full(sym)
            .ok_or(format_err!("no such element: {}", sym))?;
        let line = format!(
            "{index:5} {tag:5} {x:-18.10E} {y:-18.10E} {z:-18.10E}\n",
            index=i,
            tag = tag + 1,
            x = position[0],
            y = position[1],
            z = position[2],
        );
        lines.push_str(&line);
    }

    // format crystal
    if let Some(lat) = mol.lattice {
        let o = lat.origin();
        let line = format!("{:-18.10E} {:-18.10E} {:-18.10E}\n", o[0], o[1], o[2]);
        lines.push_str(&line);

        let vectors = lat.vectors();
        for i in 0..3 {
            let v = vectors[i];
            let line = format!("{:-18.10E} {:-18.10E} {:-18.10E}\n", v[0], v[1], v[2]);
            lines.push_str(&line);
        }
    }

    Ok(lines)
}

#[test]
fn test_dftb_molecule() {
    let mols = io::read("tests/files/mol2/LTL-crysin-ds.mol2")
        .expect("mol2 file");
    let mol = &mols[0];
    let stream = format_molecule(&mol).expect("dftb stream");
    let (_, rmol) = get_molecule(&stream).expect("dftb molecule from stream");
    assert_eq!(mol.natoms(), rmol.natoms());

    // lattice
    assert!(rmol.lattice.is_some());
    let lat1 = mol.lattice.unwrap();
    let lat2 = rmol.lattice.unwrap();
    let va1 = lat1.vector_a();
    let va2 = lat2.vector_a();
    for i in 0..3 {
        assert_relative_eq!(va1[i], va2[i], epsilon=1e-3);
    }

    // symbol and position
    for (a, b) in mol.atoms().zip(rmol.atoms()) {
        assert_eq!(a.symbol(), b.symbol());
        let pa = a.position();
        let pb = b.position();
        for i in 0..3 {
            assert_relative_eq!(pa[i], pb[i], epsilon=1e-4);
        }
    }
}
// cb010521-689c-4b07-9cba-4e2dbea7fac3 ends here

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::9f2c8e7a-461f-4f9c-aef2-636ac1f480a4][9f2c8e7a-461f-4f9c-aef2-636ac1f480a4]]
pub struct DftbInputFile();

impl ChemFileLike for DftbInputFile {
    fn ftype(&self) -> &str {
        "dftb/input"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".gen", ".dftb"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        get_molecule(&chunk)
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        format_molecule(mol)
    }
}
// 9f2c8e7a-461f-4f9c-aef2-636ac1f480a4 ends here
