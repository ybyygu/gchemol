// [[file:~/Workspace/Programming/gchemol/gchemol.note::6b59958f-7e56-4b16-b02a-cc01e5de3da8][6b59958f-7e56-4b16-b02a-cc01e5de3da8]]
use serde_json;
use serde_derive;
use indexmap::IndexMap;

use std::fs::File;

use handlebars;

use errors::*;
use molecule::Molecule;
use io;

use handlebars::{
    to_json,
    Handlebars,
    Helper,
    HelperResult,
    RenderContext,
    RenderError
};

// define a helper for formatting string or number
fn format(h: &Helper, _: &Handlebars, rc: &mut RenderContext) -> HelperResult {
    // get positional parameter from helper or throw an error
    let param = h.param(0).ok_or(RenderError::new("Param 0 is required for format helper."))?;

    // get keyword parameters
    let width =  h.hash_get("width")
        .and_then(|v| v.value().as_u64());
    let prec =  h.hash_get("prec")
        .and_then(|v| v.value().as_u64());
    let align =  h.hash_get("align")
        .and_then(|v| v.value().as_str());

    // format string
    if param.value().is_string() {
        let v = param.value()
            .as_str()
            .ok_or(RenderError::new("param 0: not str"))?;
        let width = width.unwrap_or(0) as usize;
        let rendered = if let Some(align) = align {
            match align {
                "center" => format!("{:^width$}", v, width=width),
                "right"  => format!("{:<width$}", v, width=width),
                "left"   => format!("{:>width$}", v, width=width),
                _        => format!("{:width$}", v, width=width),
            }
        } else {
            format!("{:width$}", v, width=width)
        };
        rc.writer.write(rendered.into_bytes().as_ref())?;

    // format number
    } else if param.value().is_number() || param.value().is_f64() {
        let num: f64 = param.value()
            .as_f64()
            .ok_or(RenderError::new("param 0: not f64 number"))?;

        let width = width.unwrap_or(8) as usize;
        let prec = prec.unwrap_or(4) as usize;
        let rendered = if let Some(align) = align {
            match align {
                "center" => format!("{:^width$.prec$}", num, width=width, prec=prec),
                "right"  => format!("{:<width$.prec$}", num, width=width, prec=prec),
                "left"   => format!("{:>width$.prec$}", num, width=width, prec=prec),
                _        => format!("{:width$.prec$}",  num, width=width, prec=prec),
            }
        } else {
            format!("{:-width$.prec$}", num, width=width, prec=prec)
        };
        rc.writer.write(rendered.into_bytes().as_ref())?;
    } else {
        return Err(RenderError::new("Possible type for param 0: string or number"));
    }

    Ok(())
}

#[derive(Debug, Serialize)]
struct AtomData {
    symbol: String,
    number: usize,
    x: f64,
    y: f64,
    z: f64,
    fx: f64,
    fy: f64,
    fz: f64,
}

impl Default for AtomData {
    fn default() -> Self {
        AtomData {
            symbol: "C".into(),
            number: 6,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            fx: 0.0,
            fy: 0.0,
            fz: 0.0,
        }
    }
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

#[derive(Debug, Serialize)]
struct MoleculeData {
    title: String,
    unit_cell: Option<UnitCell>,
    number_of_atoms: usize,
    number_of_bonds: usize,
    atoms: Vec<AtomData>,
    bonds: Vec<BondData>,

    // mapping element type:
    // O C H
    // 1 2 3
    element_types: IndexMap<String, usize>,
}

/// construct a shallow representation of molecule for templating
fn molecule_to_template_data(mol: &Molecule) -> serde_json::Value {
    // unit cell data
    let unit_cell = if let Some(mut lat) = mol.lattice {
        let [va, vb, vc] = lat.vectors();
        let [a, b, c] = lat.lengths();
        let [alpha, beta, gamma] = lat.angles();

        let cell = UnitCell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            va,
            vb,
            vc
        };

        Some(cell)
    } else {
        None
    };

    // atoms data
    let mut atoms = vec![];
    for a in mol.atoms() {
        let [x, y, z] = a.position();
        let number = a.number();
        let symbol = a.symbol().to_string();
        let [fx, fy, fz] = mol.lattice
            .map(|mut lat| lat.to_frac([x, y, z]))
            .unwrap_or([0.0; 3]);

        atoms.push(AtomData{
            symbol,
            number,
            x,
            y,
            z,
            fx,
            fy,
            fz
        })
    }

    let mut bonds = vec![];

    let mut element_types = indexmap!{};
    for a in mol.atoms() {
        let k = a.symbol().into();
        let c = element_types.entry(k).or_insert(0);
        *c += 1;
    }

    let mut md = MoleculeData {
        title: mol.title(),
        number_of_atoms: mol.natoms(),
        number_of_bonds: mol.nbonds(),
        unit_cell,
        atoms,
        bonds,
        element_types,
    };

    json!({
        "molecule": md,
    })
}

/// render molecule in user defined template
pub fn render_molecule_with(mol: &Molecule, template: &str) -> Result<String> {
    let data = molecule_to_template_data(mol);
    let mut h = Handlebars::new();
    h.register_helper("format", Box::new(format));
    h.render_template(template, &data).chain_err(||"failed to render")
}

#[test]
fn test_template_render() {
    let mols = io::read("tests/files/mol2/LTL-crysin-ds.mol2").expect("molecules");
    let mol = &mols[0];

    let template = io::read_file("tests/files/templates/xyz.hbs").expect("template xyz.hbs");
    let x = render_molecule_with(&mol, &template);
}
// 6b59958f-7e56-4b16-b02a-cc01e5de3da8 ends here
