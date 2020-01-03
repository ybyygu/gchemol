// use

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*use][use:1]]
use std::path::Path;
use std::fs::File;
use std::collections::HashMap;

use crate::core_utils::*;
use text_parser::*;
use nom::IResult;
use gchemol_core::{Atom, Bond, Molecule, Lattice, BondKind, AtomKind};
// use:1 ends here

// mods

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*mods][mods:1]]
mod xyz;
mod pdb;
mod mol2;
mod sdf;
mod cif;
mod vasp;
// mods:1 ends here

// traits
// Unify behaviors for all chemical file formats

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*traits][traits:1]]
pub trait ChemFileLike {
    /// file type string
    fn ftype(&self) -> &str;

    /// Supported file types in file extension, for example:
    /// [".xyz", ".mol2"]
    fn extensions(&self) -> Vec<&str>;

    /// Determine if file `filename` is parable according to its supported file extensions
    fn parsable(&self, filename: &Path) -> bool {
        let filename = format!("{}", filename.display());
        let filename = filename.to_lowercase();
        for s in self.extensions() {
            if filename.ends_with(&s.to_lowercase()) {
                return true;
            }
        }

        false
    }

    /// Formatted representation of a Molecule
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
    fn write(&self, filename: &Path, mols: &[Molecule]) -> Result<()> {
        use crate::io::prelude::ToFile;

        let txt = self.format(mols)?;
        &txt.to_file(filename)?;
        Ok(())
    }

    /// print a brief description about a chemical file format
    fn describe(&self) {
        println!("filetype: {:?}, possible extensions: {:?}",
                 self.ftype(),
                 self.extensions()
        );
    }

    /// Parse a single molecule from &str using facilities provied by nom crate
    /// file will be write-only if not implemented
    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        unimplemented!()
    }

    /// Default implementation: parse multiple molecules from `filename`
    // Note: cannot use generic type parameters
    fn parse(&self, path: &Path) -> Result<Vec<Molecule>> {
        let fp = File::open(path)
            .map_err(|e| format_err!("failed to open {}: {:?}", path.display(), e))?;

        let parser = TextParser::default();

        let mut mols = vec![];
        parser.parse(fp,
                     // parse a single part
                     |input| {
                         self.parse_molecule(input)
                     },

                     // collect all parts
                     |m| {
                         mols.push(m);
                     }
        )?;

        Ok(mols)
    }

    /// Return the last molecule in the stream
    fn parse_from(&self, stream: &str) -> Result<Molecule> {
        use ::std::io::Cursor;
        let fp = Cursor::new(stream.as_bytes());

        let mut mol = Molecule::default();
        let parser = TextParser::default();
        parser.parse(fp,
                     // parse a single part
                     |input| {
                         self.parse_molecule(input)
                     },

                     // collect all parts
                     |m| {
                         mol = m;
                     }
        )?;

        Ok(mol)
    }
}
// traits:1 ends here

// backends

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*backends][backends:1]]
macro_rules! avail_parsers {
    () => {
        vec![
            Box::new(xyz::XYZFile()),
            Box::new(xyz::PlainXYZFile()),
            Box::new(mol2::Mol2File()),
            Box::new(sdf::SdfFile()),
            Box::new(vasp::PoscarFile()),
            Box::new(cif::CifFile()),
            Box::new(pdb::PdbFile()),
        ]
    };
}

/// guess the most appropriate file format by its file extensions
pub fn guess_chemfile<P: AsRef<Path>>(path: P, fmt: Option<&str>) -> Option<Box<ChemFileLike>>{
    let filename = path.as_ref();

    let backends: Vec<Box<ChemFileLike>> = avail_parsers!();
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
pub fn describe_backends() {
    let backends: Vec<Box<ChemFileLike>> = avail_parsers!();

    for cf in backends {
        cf.describe();
    }
}
// backends:1 ends here
