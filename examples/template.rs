// [[file:~/Workspace/Programming/gchemol/gchemol.note::1189830c-e1ea-439f-acff-0bc553d89b3b][1189830c-e1ea-439f-acff-0bc553d89b3b]]
#[macro_use]
extern crate tera;
extern crate gchemol;
#[macro_use]
extern crate serde_derive;

use std::collections::HashMap;
use gchemol::errors::*;
use tera::{
    Tera,
    Context,
};

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
    let tera = Tera::new("examples/*.in").chain_err(|| "bad")?;
    let mut context = Context::new();

    let mols = gchemol::io::read("tests/files/mol2/arginyl-ds.mol2")?;
    let mol = &mols[0];
    context.add("title", &mol.title());
    context.add("natoms", &mol.natoms());
    context.add("nbonds", &mol.nbonds());
    // context.add("atoms", &vec![1,2,3,4]);
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

    context.add("atoms", &atoms);
    let rendered = tera.render("mol.in", &context).chain_err(||"Failed to render template")?;
    println!("{}", rendered);

    Ok(())
}
// 1189830c-e1ea-439f-acff-0bc553d89b3b ends here
