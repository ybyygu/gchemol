// [[file:~/Workspace/Programming/gchemol/gchemol.note::6b59958f-7e56-4b16-b02a-cc01e5de3da8][6b59958f-7e56-4b16-b02a-cc01e5de3da8]]
#[macro_use]
extern crate serde_json;
#[macro_use]
extern crate serde_derive;

extern crate gchemol;
extern crate handlebars;

use gchemol::errors::*;
use handlebars::Handlebars;

static TEMPLATE: &'static str =
"{{natoms}}
{{title}}
{{#each atoms as |a| ~}}
{{a.symbol}} {{a.x}} {{a.y}} {{a.z}}
{{/each~}}
";

#[derive(Debug, Serialize)]
struct AtomData {
    symbol: String,
    x: f64,
    y: f64,
    z: f64,
    // fx: f64,
    // fy: f64,
    // fz: f64,
}

#[derive(Debug, Serialize)]
struct MoleculeData {
    title: String,
    formula: String,
    unit_cell: UnitCell,
    natoms: usize,
    nbonds: usize,
    atoms: Vec<AtomData>,
    bonds: Vec<BondData>,
}

#[derive(Debug, Serialize)]
struct BondData {
    i: usize,
    j: usize,
    order: f64
}

#[derive(Debug, Serialize)]
struct UnitCell {
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
    va: [f64; 3],
    vb: [f64; 3],
    vc: [f64; 3],
}


fn main() -> Result<()> {
    let mut h = Handlebars::new();
    let mols = gchemol::io::read("tests/files/mol2/arginyl-ds.mol2")?;
    let mol = &mols[0];

    let mut atoms = vec![];
    for a in mol.atoms() {
        let [x, y, z] = a.position();
        let p = AtomData {
            symbol: a.symbol().to_string(),
            x,
            y,
            z
        };
        atoms.push(p);
    }

    let data = json!({
        "title": mol.title(),
        "natoms": mol.natoms(),
        "nbonds": mol.nbonds(),
        "atoms": atoms,
    });

    let x = h.render_template(TEMPLATE, &data).chain_err(|| "failed")?;
    println!("{:}", x);

    Ok(())
}
// 6b59958f-7e56-4b16-b02a-cc01e5de3da8 ends here
