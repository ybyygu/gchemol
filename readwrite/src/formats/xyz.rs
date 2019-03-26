// base

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*base][base:1]]
use super::*;
// base:1 ends here

// atom/atoms

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*atom/atoms][atom/atoms:1]]
/// Create Atom object from xyz line
/// # Example
/// C -11.4286  1.7645  0.0000
named!(pub read_atom_xyz<&str, Atom>, sp!(
    do_parse!(
        sym      : alt!(nom::alpha|nom::digit) >> // element symbol, "1" or "H"
        position : xyz_array                   >>
                   read_line                   >> // ignore the remaining characters
        (
            Atom::new(sym, position)
        )
    )
));

#[test]
fn test_parser_read_atom() {
    let (_, x) = read_atom_xyz("C -11.4286 -1.3155  0.0000\n").unwrap();
    assert_eq!("C", x.symbol());
    let (_, x) = read_atom_xyz("6 -11.4286 -1.3155  0.0000 \n").unwrap();
    assert_eq!("C", x.symbol());
    assert_eq!(0.0, x.position()[2]);
}

/// Create a list of atoms from many lines in xyz format
/// # Example
/// C -11.4286  1.7645  0.0000
/// C -10.0949  0.9945  0.0000
/// C -10.0949 -0.5455  0.0000
named!(read_atoms_xyz<&str, Vec<Atom>>, many1!(read_atom_xyz));

#[test]
fn test_parser_read_atoms() {
    let txt = "C -11.4286  1.7645  0.0000
C -10.0949  0.9945  0.0000
C -10.0949 -0.5455  0.0000
C -11.4286 -1.3155  0.0000
\n";
    let (_, atoms) = read_atoms_xyz(txt).expect("read_atoms");
    assert_eq!(4, atoms.len());
}
// atom/atoms:1 ends here

// molecule

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*molecule][molecule:1]]
/// Create a Molecule object from lines in plain xyz format (coordinates only)
named!(read_molecule_pxyz<&str, Molecule>, do_parse!(
    atoms: read_atoms_xyz >>
    (
        {
            let mut mol = Molecule::new("plain xyz");
            for a in atoms {
                mol.add_atom(a);
            }
            mol
        }
    )
));

/// Create a Molecule object from lines in xyz format
#[inline]
fn read_molecule_xyz(input: &str) -> nom::IResult<&str, Molecule> {
    // the first line
    let (input, natoms) = read_usize(input)?;
    // the second line
    let (input, title) = read_line(input)?;

    // the following lines containing coordinates
    let (input, mut mol) = read_molecule_pxyz(input)?;

    // check atoms count
    if natoms != mol.natoms() {
        warn!(
            "Malformed xyz format: expect {} atoms, but found {}",
            natoms,
            mol.natoms()
        );
    }

    // take molecule name
    if !title.trim().is_empty() {
        mol.name = title.into();
    }

    Ok((input, mol))
}

#[test]
fn test_parser_read_molecule() {
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
H -13.7062  1.5395  0.0000\n";

    let txt = format!("{}{}", txt, MAGIC_EOF);
    let (_, mol) = read_molecule_xyz(&txt).unwrap();

    assert_eq!(12, mol.natoms());
}
// molecule:1 ends here

// XYZFile

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*XYZFile][XYZFile:1]]
pub struct XYZFile();

named!(parse_xyz<&str, Molecule>, do_parse!(
    mol: read_molecule_xyz >>
    (
        mol
    )
));

impl ChemFileLike for XYZFile {
    fn ftype(&self) -> &str {
        "text/xyz"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".xyz"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        parse_xyz(chunk)
    }

    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        // meta information
        let mut lines = String::new();
        lines.push_str(&format!("{}\n", mol.natoms()));
        lines.push_str(&format!("{}\n", mol.title()));

        // coordinates
        for a in mol.atoms() {
            let p = a.position();
            let v = a.momentum();
            let sym = a.symbol();
            let s = format!(
                "{:6} {:-18.6}{:-18.6}{:-18.6}{:-18.6}{:-18.6}{:-18.6}\n",
                sym, p[0], p[1], p[2], v[0], v[1], v[2]
            );
            lines.push_str(&s);
        }

        Ok(lines)
    }
}

#[test]
fn test_formats_xyz() {
    let file = XYZFile();
    let path = Path::new("tests/files/xyz/c2h4.xyz");
    let mols = file.parse(path).expect("c2h4 xyz");
    assert_eq!(1, mols.len());
    assert_eq!(6, mols[0].natoms());

    // parse multiple molecules
    let path = Path::new("tests/files/xyz/multi.xyz");
    let mols = file.parse(path).expect("multi xyz");
    assert_eq!(6, mols.len());

    let natoms_expected = vec![16, 10, 16, 16, 16, 13];
    let natoms: Vec<_> = mols.iter().map(|m| m.natoms()).collect();
    assert_eq!(natoms_expected, natoms);
}
// XYZFile:1 ends here

// PlainXYZFile

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*PlainXYZFile][PlainXYZFile:1]]
/// plain xyz coordinates with atom symbols
#[derive(Debug, Clone)]
pub struct PlainXYZFile();

named!(parse_plain_xyz<&str, Molecule>, do_parse!(
    mol: read_molecule_pxyz >> blank_line >>
    (
        mol
    )
));

impl ChemFileLike for PlainXYZFile {
    /// possible file extensions
    fn extensions(&self) -> Vec<&str> {
        [".coord", ".pxyz", ".coords"].to_vec()
    }

    fn ftype(&self) -> &str {
        "text/pxyz"
    }

    /// Parse a single molecule from file `filename`
    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        parse_plain_xyz(chunk)
    }

    /// Return a string representation of molecule
    /// Multiple molecules will be separated by a blank line
    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        let mut lines = String::new();

        for a in mol.atoms() {
            lines.push_str(format!("{}\n", a.to_string()).as_ref());
        }

        // append a blank line as a separator between multiple molecules
        lines.push_str("\n");

        Ok(lines)
    }
}

#[test]
fn test_formats_plain_xyz() {
    let filename = Path::new("tests/files/xyz/c2h4.pxyz");
    let file = PlainXYZFile();
    assert!(file.parsable(filename));
    let mols = file.parse(filename).unwrap();
    assert_eq!(1, mols.len());
    assert_eq!(6, mols[0].natoms());

    // element numbers
    let filename = Path::new("tests/files/xyz/ele-num.pxyz");
    let mols = file.parse(filename).expect("ele-num.pyxz");
    assert_eq!(1, mols.len());
    assert_eq!(17, mols[0].natoms());
    let symbols = &mols[0].symbols();
    assert_eq!("Si", symbols[0]);

    // parse multiple molecules
    let mols = file
        .parse(Path::new("tests/files/xyz/multi.pxyz"))
        .expect("multi xyz");
    assert_eq!(6, mols.len());

    let natoms_expected = vec![16, 10, 16, 16, 16, 13];
    let natoms: Vec<_> = mols.iter().map(|m| m.natoms()).collect();
    assert_eq!(natoms_expected, natoms);
}
// PlainXYZFile:1 ends here
