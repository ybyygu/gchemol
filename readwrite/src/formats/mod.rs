// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::7faf1529-aae1-4bc5-be68-02d8ccdb9267][7faf1529-aae1-4bc5-be68-02d8ccdb9267]]
use std::str;
use std::fs::File;
use std::io::{BufRead, BufReader};

use quicli::prelude::*;
pub use gchemol_core::{
    Atom,
    Bond,
    BondKind,
    Molecule,
    Lattice,
};

/// whitespace including one or more spaces or tabs
named!(pub space_token<&str, &str>, eat_separator!(&b" \t"[..]));
macro_rules! sp (
    ($i:expr, $($args:tt)*) => (
        {
            sep!($i, space_token, $($args)*)
        }
    )
);

pub use parser::*;
pub use io;

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
pub mod dftb;

const BUF_SIZE: usize = 20 * 1024;
pub const MAGIC_EOF: &str = "$THIS_IS_THE=MAGIC_END_OF_FILE$";
// pub const MAGIC_EOF: &str = "\n\n";

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
    fn format(&self, mols: &[Molecule]) -> Result<String> {
        let mut ms = String::new();
        for mol in mols {
            let m = self.format_molecule(mol)?;
            ms.push_str(&m);
        }

        Ok(ms)
    }

    /// Save multiple molecules into a file
    fn write(&self, filename: &str, mols: &[Molecule]) -> Result<()> {
        use io::prelude::ToFile;

        let lines = self.format(mols)?;
        &lines.to_file(filename)?;
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

    // FIXME: rename to from_file
    /// Default implementation: parse multiple molecules from `filename`
    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        use nom;

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
                if chunk == MAGIC_EOF {break 'out;}

                // 1. process molecular file parsing
                match self.parse_molecule(&chunk) {
                    // 1.1 successfully parsed into one molecule
                    Ok((r, mol)) => {
                        remained.clear();
                        remained.push_str(r);
                        mols.push(mol);
                    },
                    // 1.2 when chunk is incomplete
                    // `Incomplete` means the nom parser does not have enough data to decide,
                    // so we wait for the next refill and then retry parsing
                    Err(nom::Err::Incomplete(i)) => {
                        if final_stream {
                            break 'out;
                        }
                        remained.clear();
                        remained.push_str(&chunk);
                        break;
                    },
                    // 1.3 found parse errors
                    Err(nom::Err::Error(err)) => {
                        eprintln!("found error when parsing {}: {:?}", filename, err);
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
// pub fn guess_chemfile<P: AsRef<Path>>(path: P, fmt: Option<&str>) -> Option<Box<ChemFileLike>>{
// filename: stick to &str, instead of Path
pub fn guess_chemfile(filename: &str, fmt: Option<&str>) -> Option<Box<ChemFileLike>>{
    // let path = path.as_ref();
    // let filename = &format!("{}", path.display());

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
        Box::new(dftb::DftbInputFile()),
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
            if x.parsable(filename) {
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
        Box::new(dftb::DftbInputFile()),
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