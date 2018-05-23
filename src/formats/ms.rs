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

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a2072561-bbcd-4fef-97f6-8fc67a23846a][a2072561-bbcd-4fef-97f6-8fc67a23846a]]
/// Accelrys XTL format
/// -------------------
/// The XTL file format provides a means of exchanging crystallographic
/// structure data between the Accelrys modeling products Insight II, Cerius2,
/// and Materials Studio.
pub struct XtlFile();

impl ChemFileLike for XtlFile {
    fn ftype(&self) -> &str {
        "ms/xtl"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".xtl"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        unimplemented!()
    }

    fn format(&self, mols: &Vec<Molecule>) -> Result<String> {
        if mols.len() > 1 {
            eprintln!("WARNING: only the last molecule will be kept.");
        }

        let mol = &mols.last().ok_or("no molecule")?;
        let mut lines = String::new();
        // FIXME: multiple-lines?
        lines.push_str(&format!("TITLE {}\n", mol.title()));
        lines.push_str("DIMENSION 3\n");

        if let Some(mut lattice) = mol.lattice {
            let (a, b, c) = lattice.lengths();
            let (alpha, beta, gamma) = lattice.angles();
            lines.push_str(&format!("CELL {a:.4} {b:.4} {c:.4} {alpha:2} {beta:.2} {gamma:.2}\n",
                                    a=a,
                                    b=b,
                                    c=c,
                                    alpha=alpha,
                                    beta=beta,
                                    gamma=gamma,
            ));

            lines.push_str("SYMMETRY  NUMBER 1  LABEL P1\n");
            lines.push_str("SYM MAT  1.000000  0.000000  0.000000  0.000000  1.000000  0.000000  0.000000  0.000000  1.000000 0.0000 0.0000 0.0000\n");
        } else {
            // FIXME:
            eprintln!("WARNING: no cell data found! xtl file is only suitable for solid.");
            bail!("no supported");
        }

        lines.push_str("ATOMS\n");
        lines.push_str("NAME       X          Y          Z         CHARGE   TEMP    OCCUP   SCAT\n");
        for a in mol.atoms() {
            let [x, y, z] = a.position();
            let line = format!(" {symbol:<4}{x:11.5}{y:11.5}{z:11.5}{charge:11.5}   0.0000  1.0000   {name}\n",
                               symbol=a.symbol(),
                               x=x,
                               y=y,
                               z=z,
                               name=a.symbol(),
                               charge=0.0,
            );
            lines.push_str(&line);
        }

        Ok(lines)
    }
}

#[test]
fn test_formats_xtl() {
    let mols = io::read("tests/files/mol2/LTL-crysin-ds.mol2").unwrap();
    io::write("/tmp/a.xtl", &mols);
}
// a2072561-bbcd-4fef-97f6-8fc67a23846a ends here
