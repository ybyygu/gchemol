// [[file:~/Workspace/Programming/gchemol/gchemol.note::7faf1529-aae1-4bc5-be68-02d8ccdb9267][7faf1529-aae1-4bc5-be68-02d8ccdb9267]]
use errors::*;
use io;
use Molecule;

pub trait ChemFileLike {
    /// test if file `filename` is parable
    fn parsable(&self, filename: &str) -> bool;

    /// parse molecules from file `filename`
    fn parse(&self, filename: &str) -> Result<Vec<Molecule>>;

    /// represent molecules in certain format
    fn represent(&self, mols: &Vec<Molecule>) -> String;

    fn save(&self, mols: &Vec<Molecule>, filename: &str) -> Result<()> {
        let lines = self.represent(mols);
        io::write_file(lines, filename)?;
        Ok(())
    }
}

/// plain xyz coordinates with atom symbols
pub struct PlainXYZFile();

impl ChemFileLike for PlainXYZFile {
    fn parsable(&self, filename: &str) -> bool {
        false
    }

    /// parse molecules from file `filename`
    fn parse(&self, filename: &str) -> Result<Vec<Molecule>> {
        let mut mols = vec![];

        Ok(mols)
    }

    /// Return a string representation of the last molecule in the list
    /// Return empty string if no molecule found
    fn represent(&self, mols: &Vec<Molecule>) -> String {
        let mut lines = String::new();

        if let Some(mol) = mols.last() {
            for a in mol.atoms() {
                lines.push_str(format!("{}\n", a.to_string()).as_ref());
            }
        }

        lines
    }
}

#[test]
fn test_formats_plainxyz() {
    let mol = Molecule::from_file("tests/data/c2h4.xyz").unwrap();

    let file = PlainXYZFile();
    println!("{}", file.represent(&vec![mol]));
}
// 7faf1529-aae1-4bc5-be68-02d8ccdb9267 ends here
