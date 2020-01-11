// imports

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*imports][imports:1]]
use crate::core_utils::*;
// imports:1 ends here

// header

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*header][header:1]]
//===============================================================================#
//   DESCRIPTION:  basic read & write support for molecular file
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-11 Wed 15:42>
//       UPDATED:  <2020-01-11 Sat 20:08>
//===============================================================================#
// header:1 ends here

// trait

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*trait][trait:1]]
use std::path::Path;

use crate::core_utils::*;

pub mod prelude {
    use super::*;
    use gchemol_core::Molecule;

    pub trait FromFile: Sized {
        /// Return content of text file in string
        ///
        /// don't use this if file is very large
        ///
        fn from_file<P: AsRef<Path>>(path: P) -> Result<Self>;
    }

    pub trait ToFile {
        /// write string content to an external file
        ///
        /// _Note:_ Replaces the current file content if the file already exists.
        ///
        fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<()>;
    }

    pub trait StringIO {
        /// format molecule as string in specific type
        fn format_as<S: AsRef<str>>(&self, fmt: S) -> Result<String>;

        /// construct molecule from string in specific type
        fn parse_from<S: AsRef<str>, T: AsRef<str>>(s: S, fmt: T) -> Result<Molecule>;
    }
}

// FIXME: to be removed.
// reexport
pub use guts::fs::read_file;

impl FromFile for String {
    fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        guts::fs::read_file(path)
    }
}

impl ToFile for str {
    fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        guts::fs::write_to_file(path, &self)
    }
}
// trait:1 ends here

// molecule

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*molecule][molecule:1]]
use gchemol_core::{Atom, Molecule};

// import important traits
use crate::io::prelude::*;

impl FromFile for Molecule {
    /// Construct molecule from external text file
    fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let cf = guess_chemfile_from_filename(path)?;
        let mut mols = cf.parse(path)?;
        mols.pop()
            .ok_or(format_err!("No molecule: {:?}", path.display()))
    }
}

impl ToFile for Molecule {
    /// Save molecule to an external file
    fn to_file<T: AsRef<Path>>(&self, filename: T) -> Result<()> {
        let filename = filename.as_ref();
        let cf = guess_chemfile(&filename.display().to_string(), None)
            .ok_or(format_err!("not supported file format: {:?}", filename))?;
        let t = cf.format_molecule(&self)?;

        t.to_file(filename)?;

        Ok(())
    }
}

impl StringIO for Molecule {
    /// format molecule as string in specific type
    fn format_as<S: AsRef<str>>(&self, fmt: S) -> Result<String> {
        let fmt = fmt.as_ref();
        let cf = guess_chemfile_from_fmt(fmt)?;
        cf.format_molecule(&self)
    }

    // fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
    /// construct molecule from string in specific type
    fn parse_from<S: AsRef<str>, T: AsRef<str>>(s: S, fmt: T) -> Result<Molecule> {
        let fmt = fmt.as_ref();
        let cf = guess_chemfile_from_fmt(fmt)?;

        let m = cf.parse_from(s.as_ref())?;

        Ok(m)
    }
}

fn guess_chemfile_from_fmt(fmt: &str) -> Result<Box<ChemFileLike>> {
    let msg = format_err!("not supported format: {:?}", fmt);
    guess_chemfile("", Some(fmt)).ok_or(msg)
}

fn guess_chemfile_from_filename<P: AsRef<Path>>(path: P) -> Result<Box<ChemFileLike>> {
    let filename = format!("{}", path.as_ref().display());
    let msg = format_err!("not supported format: {:?}", filename);
    guess_chemfile(&filename, None).ok_or(msg)
}

#[test]
fn test_molecule_formats() {
    let mol = Molecule::from_file("tests/files/xyz/c2h4.xyz").unwrap();
    assert_eq!(6, mol.natoms());

    let txt = mol.format_as("text/xyz").expect("format mol as xyz");
    let mol2 = Molecule::parse_from(txt, "text/xyz").expect("parse molecule from xyz stream");
    assert_eq!(mol.natoms(), mol2.natoms());
}
// molecule:1 ends here

// molecules

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*molecules][molecules:1]]
use crate::formats::{guess_chemfile, ChemFileLike};

type ReadResult = Result<Vec<Molecule>>;

/// Some options for reading/writing chemical file
pub struct FileOptions<'a> {
    chemfile: Option<Box<ChemFileLike>>,
    fmt: Option<&'a str>,
}

// impl<T: ChemFileLike> Default for FileOptions<T> {
impl<'a> Default for FileOptions<'a> {
    fn default() -> Self {
        FileOptions {
            chemfile: None,
            fmt: None,
        }
    }
}

// impl<T: ChemFileLike> FileOptions<T> {
impl<'a> FileOptions<'a> {
    pub fn new() -> Self {
        FileOptions::default()
    }

    /// set chemfile format
    pub fn fmt(&mut self, fmt: &'a str) -> &mut Self {
        self.fmt = Some(fmt);
        self
    }

    /// read molecules from file
    pub fn read<P: AsRef<Path>>(&self, path: P) -> ReadResult {
        let path = path.as_ref();
        let msg = format_err!(
            "not supported file\nfilename: {:}\n fmt: {:?}",
            path.display(),
            self.fmt
        );
        let chemfile = guess_chemfile(path, self.fmt).ok_or(msg)?;
        let mols = chemfile.parse(path)?;

        Ok(mols)
    }

    pub fn write<P: AsRef<Path>>(&self, path: P, mols: &[Molecule]) -> Result<()> {
        let path = path.as_ref();
        let msg = format_err!(
            "not supported file\nfilename: {:}\n fmt: {:?}",
            path.display(),
            self.fmt
        );
        let chemfile = guess_chemfile(path, self.fmt).ok_or(msg)?;
        chemfile.write(path, mols)
    }
}

/// Read an iterator over `Molecule` from file.
/// file format will be determined according to the path
// https://stackoverflow.com/questions/26368288/how-do-i-stop-iteration-and-return-an-error-when-iteratormap-returns-a-result
pub fn read<P: AsRef<Path>>(path: P) -> Result<Vec<Molecule>> {
    let path = path.as_ref();
    FileOptions::new().read(path)
}

/// Write molecules into file
/// file format will be determined according to the path
pub fn write<P: AsRef<Path>>(path: P, mols: &[Molecule]) -> Result<()> {
    let path = path.as_ref();
    FileOptions::new().write(path, mols)
}

#[test]
fn test_io_read_plain_xyz() {
    let mols = FileOptions::new()
        .fmt("text/pxyz")
        .read("tests/files/xyz/ele-num.pxyz")
        .unwrap();

    assert_eq!(1, mols.len());
    assert_eq!(17, mols[0].natoms());
}
// molecules:1 ends here

// trajectory

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*trajectory][trajectory:1]]
use crate::gchemol_core::trajectory::*;
use std::convert::TryFrom;
use std::convert::TryInto;

impl FromFile for Trajectory {
    /// Construct molecule from external text file
    fn from_file<P: AsRef<Path>>(path: P) -> Result<Trajectory> {
        let path = path.as_ref();
        guess_chemfile_from_filename(path)?.parse(path)?.try_into()
    }
}

impl ToFile for Trajectory {
    /// Save trajectory to an external file
    fn to_file<T: AsRef<Path>>(&self, filename: T) -> Result<()> {
        let filename = filename.as_ref();

        let mols: Vec<_> = self.iter().collect();
        write(filename, &mols)
    }
}
// trajectory:1 ends here
