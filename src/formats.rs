// [[file:~/Workspace/Programming/gchemol/gchemol.note::7faf1529-aae1-4bc5-be68-02d8ccdb9267][7faf1529-aae1-4bc5-be68-02d8ccdb9267]]
use errors::*;
use io;
use std::path::{Path, PathBuf};
use Atom;
use Molecule;

pub trait ChemFileLike {
    /// file type string
    fn ftype(&self) -> &str;

    /// Supported file types in file extension, for example:
    /// [".xyz", ".mol2"]
    fn extensions(&self) -> Vec<&str>;

    /// test if file `filename` is parable
    fn parsable<P: AsRef<Path>>(&self, path: P) -> bool {
        let path = format!("{}", path.as_ref().display());
        let path = path.to_lowercase();
        for s in self.extensions() {
        if path.ends_with(&s.to_lowercase()) {
                return true;
            }
        }

        false
    }

    /// parse molecules from file `filename`
    fn parse<P: AsRef<Path>>(&self, filename: P) -> Result<Vec<Molecule>>;

    /// format molecules in certain format
    fn format(&self, mols: &Vec<Molecule>) -> Result<String> {
        bail!("No implemented yet.");
    }

    /// Save multiple molecules in a file
    fn save(&self, mols: &Vec<Molecule>, filename: &str) -> Result<()> {
        let lines = self.format(mols)?;
        io::write_file(lines, filename)?;
        Ok(())
    }
}

/// plain xyz coordinates with atom symbols
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
    fn parse<P: AsRef<Path>>(&self, filename: P) -> Result<Vec<Molecule>> {
        let path = filename.as_ref();
        let txt = io::read_file(path)?;
        let mut mol = Molecule::new("from plain coordinates");
        for line in txt.lines() {
            let line = line.trim();
            if line.is_empty() {
                break;
            }
            let a: Atom = line.parse()?;
            mol.add_atom(a);
        }

        Ok([mol].to_vec())
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
// 7faf1529-aae1-4bc5-be68-02d8ccdb9267 ends here
