// [[file:~/Workspace/Programming/gchemol/gchemol.note::7faf1529-aae1-4bc5-be68-02d8ccdb9267][7faf1529-aae1-4bc5-be68-02d8ccdb9267]]
use io;
use std::str;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub use parser::*;
use quicli::prelude::*;

pub use Atom;
pub use Molecule;
pub use lattice::Lattice;
pub use Bond;
pub use BondKind;

pub mod template;
pub mod xyz;
pub mod mol2;
pub mod pdb;
pub mod vasp;
pub mod sdf;
pub mod cif;
pub mod gaussian;
pub mod ms;
pub mod siesta;
pub mod gulp;

const BUF_SIZE: usize = 20 * 1024;
pub const MAGIC_EOF: &str = "$THIS_IS_THE=MAGIC_END_OF_FILE$";

use nom;

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
        bail!("unimplemented yet");
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

    /// Save multiple molecules into a file
    fn write(&self, filename: &str, mols: &Vec<Molecule>) -> Result<()> {
        let lines = self.format(mols)?;
        io::write_file(lines, filename)?;
        Ok(())
    }

    /// print a brief description about a chemical file format
    fn describe(&self) {
        println!("filetype: {:?}, possible extensions: {:?}",
                 self.ftype(),
                 self.extensions());
    }

    /// Parse a single molecule from &str using facilities provied by nom crate
    /// file will be write-only if not implemented
    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        unimplemented!()
    }

    /// Default implementation: parse multiple molecules from `filename`
    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        let fp = File::open(filename)
            .map_err(|e| format_err!("failed to open {}: {:?}", filename, e))?;

        let mut reader = BufReader::with_capacity(BUF_SIZE, fp);

        let mut mols: Vec<Molecule> = vec![];
        let mut remained = String::new();
        let mut chunk = String::new();
        let mut final_stream = false;

        // streaming the file parsing
        // FIXME: how to parse a binary file?
        // FIXME: if the file is extremely large
        'out: loop {
            // restrict reader/buffer variables scope
            let (new, length) = {
                let buffer = reader.fill_buf()?;
                let length = buffer.len();

                // a workaround for nom 4.0 changes: append a magic eof to make stream `complete`
                let new: String = if length == 0 {
                    final_stream = true;
                    String::from(MAGIC_EOF)
                } else {
                    str::from_utf8(&buffer)?.to_string()
                };

                (
                    new,
                    length,
                )
            };

            // 0. fill the chunk
            // chunk = remained + buffer
            chunk.clear();
            // remained data by nom parser
            chunk.push_str(&remained);
            // new data from the file

            // FIXME: define a maximal size to which the buffer can grow to
            // avoid out of memory death
            chunk.push_str(&new);

            loop {
                // 1. process molecular file parsing
                match self.parse_molecule(&chunk) {
                    // 1.1 successfully parsed into one molecule
                    Ok((r, mol)) => {
                        // println!("ooooooooooooooooo {}", filename);
                        // println!("{:#}", &chunk);
                        // println!("--ooooooooooooooo {}", filename);
                        // println!("{:#}", r);
                        remained.clear();
                        remained.push_str(r);
                        mols.push(mol);
                    },
                    // 1.2 when chunk is incomplete
                    // `Incomplete` means the nom parser does not have enough data to decide,
                    // so we wait for the next refill and then retry parsing
                    Err(nom::Err::Incomplete(i)) => {
                        // println!("xxxxxxxxinnnnnnnnnnnnn {}", filename);
                        // println!("{:#}", &chunk);
                        // // should not happen
                        // if final_stream {
                        //     println!("nom parser warning: fixmefixmefixme");
                        // }
                        remained.clear();
                        remained.push_str(&chunk);
                        break;
                    },
                    // 1.3 found parse errors
                    Err(nom::Err::Error(err)) => {
                        eprintln!("found error in {}: {:?}", filename, err);
                        break 'out;
                    },
                    // 1.4 found serious errors
                    Err(nom::Err::Failure(err)) => {
                        bail!("hard failure in {}: {:?}", filename, err);
                    },
                }
                // clear chunk
                chunk.clear();
                chunk.push_str(&remained);
            }

            if final_stream {
                break
            }

            // consume the reader buffer
            reader.consume(length);
        }

        Ok(mols)
    }
}


/// guess the most appropriate file format by its file extensions
pub fn guess_chemfile(path: &str, fmt: Option<&str>) -> Option<Box<ChemFileLike>>{
    let backends: Vec<Box<ChemFileLike>> = vec![
        Box::new(xyz::XYZFile()),
        Box::new(xyz::PlainXYZFile()),
        Box::new(mol2::Mol2File()),
        Box::new(sdf::SdfFile()),
        Box::new(vasp::PoscarFile()),
        Box::new(cif::CifFile()),
        Box::new(pdb::PdbFile()),
        Box::new(gaussian::GaussInputFile()),
        Box::new(ms::CarFile()),
        Box::new(ms::XtlFile()),
        Box::new(siesta::FdfFile()),
        Box::new(gulp::GulpInputFile()),
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
        Box::new(xyz::PlainXYZFile()),
        Box::new(mol2::Mol2File()),
        Box::new(sdf::SdfFile()),
        Box::new(vasp::PoscarFile()),
        Box::new(cif::CifFile()),
        Box::new(pdb::PdbFile()),
        Box::new(gaussian::GaussInputFile()),
        Box::new(ms::CarFile()),
        Box::new(ms::XtlFile()),
        Box::new(siesta::FdfFile()),
        Box::new(gulp::GulpInputFile()),
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
