// [[file:~/Workspace/Programming/gchemol/gchemol.note::a05dc011-24db-4b5c-8ae8-502d39007444][a05dc011-24db-4b5c-8ae8-502d39007444]]
use super::*;

/// Accelrys/MSI Biosym/Insight II CAR format
pub struct CarFile();

impl ChemFileLike for CarFile {
    fn ftype(&self) -> &str {
        "ms/car"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".car"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        unimplemented!()
    }

    fn format(&self, mols: &Vec<Molecule>) -> Result<String> {
        if mols.len() > 1 {
            eprintln!("WARNING: only the last molecule will be used.");
        }

        let mol = &mols.last().ok_or("no molecule")?;
        let mut lines = String::new();
        lines.push_str("!BIOSYM archive 3\nPBC=OFF\n\n");
        lines.push_str("!xxxx\n");

        for a in mol.atoms() {
            let [x, y, z] = a.position();
            let line = format!("{symbol:5}{x:15.9}{y:15.9}{z:15.9}{spc:21}{name:4}{charge:-5.3}\n",
                               symbol=a.symbol(),
                               x=x,
                               y=y,
                               z=z,
                               spc=" ",
                               name=a.symbol(),
                               charge=0.0,
            );
            lines.push_str(&line);
        }
        lines.push_str("end\nend\n");

        Ok(lines)
    }
}

#[test]
#[ignore]
fn test_formats_car() {
    let mols = io::read("tests/files/mol2/arginyl-ds.mol2").unwrap();
    io::write("/tmp/a.car", &mols);
}
// a05dc011-24db-4b5c-8ae8-502d39007444 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::78da2e68-530b-4b61-9fc8-6b0b08c2bc6e][78da2e68-530b-4b61-9fc8-6b0b08c2bc6e]]
/// Accelrys/XSD input format
pub struct XsdFile();
// 78da2e68-530b-4b61-9fc8-6b0b08c2bc6e ends here
