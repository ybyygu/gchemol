// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::cb1bceac-b835-44db-b09b-2944a6a31b32][cb1bceac-b835-44db-b09b-2944a6a31b32]]
use super::*;
use indexmap;

/// represent a siesta FDF input file(write only)
pub struct FdfFile();

impl ChemFileLike for FdfFile {
    fn ftype(&self) -> &str {
        "siesta/fdf"
    }

    fn extensions(&self) -> Vec<&str> {
        vec![".fdf"]
    }

    fn parse_molecule<'a>(&self, chunk: &'a str) -> IResult<&'a str, Molecule> {
        unimplemented!()
    }

    fn format(&self, mols: &[Molecule]) -> Result<String> {
        if mols.len() > 1 {
            eprintln!("WARNING: only the last molecule will be used.");
        }

        let mol = &mols.last().ok_or(format_err!("no molecule"))?;
        let mut lines = String::new();
        lines.push_str(&format!("SystemLabel               {}\n", "siesta"));
        lines.push_str(&format!("NumberOfAtoms             {}\n", mol.natoms()));

        // define elemental tags
        // let mut specs = std::collections::HashSet::new();
        let mut specs = indexmap::IndexSet::new();
        for a in mol.atoms() {
            specs.insert((a.number(), a.symbol()));
        }
        lines.push_str(&format!("NumberOfSpecies           {}\n", specs.len()));
        lines.push_str("%block ChemicalSpeciesLabel\n");

        let mut i = 1;
        for (n, s) in specs.iter() {
            lines.push_str(&format!(" {index:<3}{number:4}   {symbol:5>}\n",
                                    index=i,
                                    number=n,
                                    symbol=s));
            i += 1;
        }
        lines.push_str("%endblock ChemicalSpeciesLabel\n");
        lines.push_str("AtomicCoordinatesFormat   NotScaledCartesianAng\n");
        lines.push_str("%block AtomicCoordinatesAndAtomicSpecies\n");

        for a in mol.atoms() {
            let (i, _) = specs.get_full(&(a.number(), a.symbol())).expect("fdf (number, symbol) index");
            let tag = 1 + i;
            let [x, y, z] = a.position();
            let line = format!(" {x:-12.6}{y:-12.6}{z:-12.6}  {tag}\n",
                               x=x,
                               y=y,
                               z=z,
                               tag=tag,
            );
            lines.push_str(&line);
        }
        lines.push_str("%endblock AtomicCoordinatesAndAtomicSpecies\n");
        lines.push_str("LatticeConstant            1.0 ang\n");
        lines.push_str("%block LatticeVectors\n");

        if let Some(lattice) = mol.lattice {
            for [x, y, z] in lattice.vectors().into_iter() {
                lines.push_str(&format!("{:7.3}{:7.3}{:7.3}\n", x, y, z));
            }
        } else {
            lines.push_str("# !!please define lattice vectors here !!\n");
        }
        lines.push_str("%endblock LatticeVectors\n");
        lines.push_str("# reading more FDF information from default.fdf\n");
        lines.push_str("%include default.fdf\n");

        // final blank lines
        lines.push_str("\n");

        Ok(lines)
    }
}

#[test]
#[ignore]
fn test_formats_fdf() {
    let mols = io::read("tests/files/mol2/arginyl-ds.mol2").unwrap();
    io::write("/tmp/a.fdf", &mols).unwrap();
}
// cb1bceac-b835-44db-b09b-2944a6a31b32 ends here
