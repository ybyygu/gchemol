// [[file:~/Workspace/Programming/gchemol/gchemol.note::7faf1529-aae1-4bc5-be68-02d8ccdb9267][7faf1529-aae1-4bc5-be68-02d8ccdb9267]]
use errors::*;
use io;
use std::str;
use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufRead, BufReader};

pub use nom::IResult;
pub use Atom;
pub use Molecule;
pub use lattice::Lattice;

pub use parser::{
    space,
    space_token,
    take_until_end_of_line,
    digit,
    double_s,
    alphanumeric,
    xyz_array,
    alpha,
    signed_digit,
    unsigned_digit,
    end_of_line,
    not_line_ending,
    line_ending,
};

pub mod xyz;
pub mod mol2;
pub mod pdb;
pub mod vasp;
pub mod sdf;
pub mod cif;
pub mod gaussian;

const BUF_SIZE: usize = 8 * 1024;

/// Unified behaviors for all chemical file formats
pub trait ChemFileLike {
    /// file type string
    fn ftype(&self) -> &str;

    /// Supported file types in file extension, for example:
    /// [".xyz", ".mol2"]
    fn extensions(&self) -> Vec<&str>;

    /// Determine if file `filename` is parable according to its supported file extensions
    fn parsable(&self, filename: &str) -> bool {
        let filename = filename.to_lowercase();
        for s in self.extensions() {
            if filename.ends_with(&s.to_lowercase()) {
                return true;
            }
        }

        false
    }

    /// Formatted representation of a molecule
    fn format_molecule(&self, mol: &Molecule) -> Result<String> {
        unimplemented!()
    }

    /// format molecules in certain format
    /// file will be read-only if not implemented
    fn format(&self, mols: &Vec<Molecule>) -> Result<String> {
        let mut ms = String::new();
        for mol in mols {
            let m = self.format_molecule(mol)?;
            ms.push_str(&m);
        }

        Ok(ms)
    }

    /// Save multiple molecules in a file
    fn to_file(&self, mols: &Vec<Molecule>, filename: &str) -> Result<()> {
        let lines = self.format(mols)?;
        io::write_file(lines, filename)?;
        Ok(())
    }

    /// brief description about a chemical file format
    fn describe(&self) {
        println!("filetype: {:?}, possible extensions: {:?}",
                 self.ftype(),
                 self.extensions());
    }

    /// Parse a single molecule from &str using facilities provied by nom crate
    /// file will be write-only if not implemented
    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        unimplemented!();
    }

    /// Default implementation: parse multiple molecules from `filename`
    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        let fp = File::open(filename).chain_err(|| "failed")?;
        let mut reader = BufReader::with_capacity(BUF_SIZE, fp);

        let mut mols: Vec<Molecule> = vec![];
        let mut remained = String::new();
        let mut chunk = String::new();
        'out: loop {
            let length = {
                let buffer = reader.fill_buf().chain_err(|| "file buffer reading error")?;
                let new = str::from_utf8(&buffer).unwrap().to_string();
                // chunk = buffer + remained
                chunk.clear();
                chunk.push_str(&remained);
                chunk.push_str(&new);
                loop {
                    // println!("inner {:?}", j);
                    // fill chunk with remained data
                    match self.parse_molecule(&chunk) {
                        IResult::Error(err) => {
                            eprintln!("{:?}", err);
                            eprintln!("{:}", chunk);
                            break 'out;
                        },
                        IResult::Done(r, mol) => {
                            // println!("got mol with {:?} atoms", mol.natoms());
                            mols.push(mol);
                            remained = String::from(r);
                        },
                        IResult::Incomplete(i) => {
                            // eprintln!("need data: {:?}", i);
                            remained = chunk.clone();
                            break
                        },
                    }
                    // clear chunk
                    chunk.clear();
                    chunk.push_str(&remained);
                }
                buffer.len()
            };

            if length == 0 { break; }
            reader.consume(length);
        }

        Ok(mols)
    }
}

/// plain xyz coordinates with atom symbols
#[derive(Debug)]
pub struct PlainXYZFile();

impl ChemFileLike for PlainXYZFile {
    /// possible file extensions
    fn extensions(&self) -> Vec<&str> {
        [".coord"].to_vec()
    }

    fn ftype(&self) -> &str {
        "text/coord"
    }

    /// parse molecules from file `filename`
    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        let txt = io::read_file(filename)?;
        let mut mol = Molecule::new("from plain coordinates");
        for line in txt.lines() {
            let line = line.trim();
            if line.is_empty() {
                break;
            }
            let a: Atom = line.parse()?;
            mol.add_atom(a);
        }

        Ok(vec![mol])
    }

    /// Return a string representation of the last molecule in the list
    /// Return empty string if no molecule found
    fn format(&self, mols: &Vec<Molecule>) -> Result<String> {
        let mut lines = String::new();

        if let Some(mol) = mols.last() {
            for a in mol.atoms() {
                lines.push_str(format!("{}\n", a.to_string()).as_ref());
            }
        }

        Ok(lines)
    }
}

#[test]
fn test_formats_plainxyz() {
    let filename = "tests/files/plain-coords/test.coord";
    let file = PlainXYZFile();
    assert!(file.parsable(filename));
    let mols = file.parse(filename).unwrap();
    assert_eq!(1, mols.len());
    assert_eq!(12, mols[0].natoms());
}

/// guess the most appropriate file format by its file extensions
pub fn guess_chemfile(path: &str, fmt: Option<&str>) -> Option<Box<ChemFileLike>>{
    let backends: Vec<Box<ChemFileLike>> = vec![
        Box::new(xyz::XYZFile()),
        Box::new(PlainXYZFile()),
        Box::new(mol2::Mol2File()),
        Box::new(vasp::POSCARFile()),
        Box::new(cif::CifFile()),
    ];

    // 1. by file type
    if let Some(fmt) = fmt {
        for x in backends {
            if x.ftype() == fmt.to_lowercase() {
                return Some(x);
            }
        }
    // 2. or by file extension
    } else {
        for x in backends {
            if x.parsable(path) {
                return Some(x);
            }
        }
    }

    // 3. return None if no suitable backend
    None
}

/// description of all backends
// FIXME:
pub fn describe_backends() {
    let backends: Vec<Box<ChemFileLike>> = vec![
        Box::new(xyz::XYZFile()),
        Box::new(mol2::Mol2File()),
        Box::new(PlainXYZFile()),
        Box::new(vasp::POSCARFile()),
        Box::new(cif::CifFile()),
    ];

    for cf in backends {
        cf.describe();
    }
}

#[test]
#[ignore]
fn test_formats_descrb() {
    describe_backends();
}
// 7faf1529-aae1-4bc5-be68-02d8ccdb9267 ends here
