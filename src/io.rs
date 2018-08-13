// [[file:~/Workspace/Programming/gchemol/gchemol.note::891f59cf-3963-4dbe-a7d2-48279723b72e][891f59cf-3963-4dbe-a7d2-48279723b72e]]
//===============================================================================#
//   DESCRIPTION:  basic read & write support for molecular file
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-11 Wed 15:42>
//       UPDATED:  <2018-08-13 Mon 11:55>
//===============================================================================#
// 891f59cf-3963-4dbe-a7d2-48279723b72e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::4a0ea759-222c-42b7-ab82-8eb969751781][4a0ea759-222c-42b7-ab82-8eb969751781]]
use std::path::Path;

use quicli;
use quicli::prelude::*;

pub mod prelude {
    use super::*;

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
}

// reexport
pub use quicli::fs::read_file;

impl FromFile for String {
    fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        quicli::fs::read_file(path)
    }
}

impl ToFile for str {
    fn to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        quicli::fs::write_to_file(path, &self)
    }
}
// 4a0ea759-222c-42b7-ab82-8eb969751781 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::00093c10-2247-4242-a287-e5640c00cadb][00093c10-2247-4242-a287-e5640c00cadb]]
use {
    Atom,
    Molecule,
};

// import important traits
use io::prelude::*;

// fn file_extension_lower(path: &Path) -> Result<String> {
//     let ext = path.extension().ok_or(format_err!("cannot find file extension"))?;
//     let ext = ext.to_str().ok_or(format_err!("cannot handle wield file extention"))?;
//     let ext = ext.to_lowercase();

//     Ok(ext.to_string())
// }

impl FromFile for Molecule {
    /// Construct molecule from external text file
    fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let cf = guess_chemfile_from_filename(path)?;
        let filename = format!("{}", path.display());
        // FIXME: Path or &str?
        let mut mols = cf.parse(&filename)?;
        mols.pop().ok_or(format_err!("no molecule: {:?}", path.display()))
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

impl Molecule {
    /// format molecule as string in specific type
    pub fn format_as<S: AsRef<str>>(&self, fmt: S) -> Result<String> {
        let fmt = fmt.as_ref();
        let cf = guess_chemfile_from_fmt(fmt)?;
        cf.format_molecule(&self)
    }

    // fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
    /// construct molecule from string in specific type
    pub fn parse_from<S: AsRef<str>, T: AsRef<str>>(s: S, fmt: T) -> Result<Molecule> {
        let fmt = fmt.as_ref();
        let cf = guess_chemfile_from_fmt(fmt)?;

        let s = s.as_ref();
        // FIXME: ....
        let (_, m) = cf.parse_molecule(s).map_err(|e| format_err!("{:?}", e))?;

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
// 00093c10-2247-4242-a287-e5640c00cadb ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::d317857c-a18d-4630-9155-119cf3533d5a][d317857c-a18d-4630-9155-119cf3533d5a]]
use formats::{
    ChemFileLike,
    guess_chemfile,
};

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
        self.fmt = Some(fmt); self
    }

    /// read molecules from file
    pub fn read<P: AsRef<Path>>(&self, path: P) -> ReadResult {
        let path = path.as_ref();
        let filename = &format!("{}", path.display());
        let msg = format_err!("not supported file\nfilename: {:}\n fmt: {:?}", filename, self.fmt);
        let chemfile = guess_chemfile(filename, self.fmt).ok_or(msg)?;
        let mols = chemfile.parse(filename)?;

        Ok(mols)
    }

    pub fn write<P: AsRef<Path>>(&self, path: P, mols: &Vec<Molecule>) -> Result<()> {
        let path = path.as_ref();
        let filename = &format!("{}", path.display());
        let msg = format_err!("not supported file\nfilename: {:}\n fmt: {:?}", filename, self.fmt);
        let chemfile = guess_chemfile(filename, self.fmt).ok_or(msg)?;
        chemfile.write(filename, mols)
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
pub fn write<P: AsRef<Path>>(path: P, mols: &Vec<Molecule>) -> Result<()>{
    let path = path.as_ref();
    FileOptions::new().write(path, mols)
}

#[test]
fn test_io_read_coord() {
    let mols = FileOptions::new()
        .fmt("text/coord")
        .read("tests/files/plain-coords/test.xyz").unwrap();

    assert_eq!(1, mols.len());
    assert_eq!(17, mols[0].natoms());
}
// d317857c-a18d-4630-9155-119cf3533d5a ends here
