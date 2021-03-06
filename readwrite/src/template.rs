// traits

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*traits][traits:1]]
/// Render molecule in user defined format
pub trait TemplateRendering {
    /// Render molecule with user defined template
    fn render_with(&self, template: &str) -> Result<String>;
}

impl TemplateRendering for Molecule {
    fn render_with(&self, template: &str) -> Result<String> {
        render_molecule_with(&self, &template)
    }
}
// traits:1 ends here

// imports

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*imports][imports:1]]
use handlebars::*;
use serde_derive;

use indexmap::{indexmap, IndexMap};
use std::fs::File;

use crate::core_utils::*;
use crate::io;

use gchemol_core::{Atom, Molecule};
// imports:1 ends here

// format float number

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*format%20float%20number][format float number:1]]
// https://docs.rs/handlebars/1.0.0/handlebars/trait.HelperDef.html
// define a helper for formatting string or number
fn format(
    h: &Helper,
    _: &Handlebars,
    _: &Context,
    rc: &mut RenderContext,
    out: &mut dyn Output,
) -> HelperResult {
    // get positional parameter from helper or throw an error
    let param = h
        .param(0)
        .ok_or(RenderError::new("Param 0 is required for format helper."))?;

    // get keyword parameters
    let width = h.hash_get("width").and_then(|v| v.value().as_u64());
    let prec = h.hash_get("prec").and_then(|v| v.value().as_u64());
    let align = h.hash_get("align").and_then(|v| v.value().as_str());

    // format string
    if param.value().is_string() {
        let v = param
            .value()
            .as_str()
            .ok_or(RenderError::new("param 0: not str"))?;
        let width = width.unwrap_or(0) as usize;
        let rendered = if let Some(align) = align {
            match align {
                "center" => format!("{:^width$}", v, width = width),
                "right" => format!("{:<width$}", v, width = width),
                "left" => format!("{:>width$}", v, width = width),
                _ => format!("{:width$}", v, width = width),
            }
        } else {
            format!("{:width$}", v, width = width)
        };
        out.write(rendered.as_ref())?;

    // format number
    } else if param.value().is_number() || param.value().is_f64() {
        let num: f64 = param
            .value()
            .as_f64()
            .ok_or(RenderError::new("param 0: not f64 number"))?;

        let width = width.unwrap_or(18) as usize;
        let prec = prec.unwrap_or(8) as usize;
        let rendered = if let Some(align) = align {
            match align {
                "center" => format!("{:^width$.prec$}", num, width = width, prec = prec),
                "right" => format!("{:<width$.prec$}", num, width = width, prec = prec),
                "left" => format!("{:>width$.prec$}", num, width = width, prec = prec),
                _ => format!("{:width$.prec$}", num, width = width, prec = prec),
            }
        } else {
            format!("{:-width$.prec$}", num, width = width, prec = prec)
        };
        out.write(rendered.as_ref())?;
    } else {
        return Err(RenderError::new(
            "Possible type for param 0: string or number",
        ));
    }

    Ok(())
}
// format float number:1 ends here

// core

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*core][core:1]]
#[derive(Debug, Serialize)]
struct AtomData {
    index: usize,
    element_index: usize,
    symbol: String,
    number: usize,
    x: f64,
    y: f64,
    z: f64,
    fx: f64,
    fy: f64,
    fz: f64,
    vx: f64,
    vy: f64,
    vz: f64,
}

impl Default for AtomData {
    fn default() -> Self {
        AtomData {
            index: 0,
            element_index: 0,
            symbol: "C".into(),
            number: 6,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            fx: 0.0,
            fy: 0.0,
            fz: 0.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
        }
    }
}

#[derive(Debug, Serialize)]
struct BondData {
    i: usize,
    j: usize,
    order: f64,
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
struct SpeciesData {
    index: usize,
    element_symbol: String,
    element_number: usize,
    number_of_atoms: usize,
}

#[derive(Debug, Serialize)]
struct MoleculeData {
    title: String,
    unit_cell: Option<UnitCell>,
    number_of_atoms: usize,
    number_of_bonds: usize,
    number_of_species: usize,
    atoms: Vec<AtomData>,
    bonds: Vec<BondData>,

    // mapping element type:
    // O C H
    // 1 2 3
    element_types: Vec<(String, usize)>,
    species: Vec<SpeciesData>,
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
            vc,
        };

        Some(cell)
    } else {
        None
    };

    let mut bonds = vec![];

    let mut element_types: IndexMap<String, usize> = indexmap!{};
    for a in mol.atoms() {
        let k = a.symbol().into();
        let c = element_types.entry(k).or_insert(0);
        *c += 1;
    }

    // atoms data
    let mut atoms = vec![];
    for (i, a) in mol.view_atoms() {
        let [x, y, z] = a.position();
        let index = i;
        let number = a.number();
        let symbol = a.symbol().to_string();
        let [fx, fy, fz] = mol
            .lattice
            .map(|mut lat| lat.to_frac([x, y, z]))
            .unwrap_or([0.0; 3]);
        let element_index = {
            let (x, _, _) = element_types
                .get_full(a.symbol())
                .expect("element type index");
            x + 1
        };

        let [vx, vy, vz] = a.momentum();
        atoms.push(AtomData {
            index,
            element_index,
            symbol,
            number,
            x,
            y,
            z,
            fx,
            fy,
            fz,
            vx,
            vy,
            vz,
        })
    }

    // convert indexmap to plain list
    let element_types: Vec<(_, _)> = element_types.into_iter().collect();

    let n = element_types.len();
    let species: Vec<_> = element_types
        .iter()
        .enumerate()
        .map(|(i, (s, n))| SpeciesData {
            index: i + 1,
            element_symbol: s.clone(),
            // FIXME: dirty
            element_number: Atom::new(s, [0.0; 3]).number(),
            number_of_atoms: *n,
        })
        .collect();

    let md = MoleculeData {
        title: mol.title(),
        number_of_atoms: mol.natoms(),
        number_of_bonds: mol.nbonds(),
        number_of_species: n,
        unit_cell,
        atoms,
        bonds,
        element_types,
        species,
    };

    json!({
        "molecule": md,
    })
}

/// render molecule in user defined template
pub fn render_molecule_with(mol: &Molecule, template: &str) -> Result<String> {
    let mut h = Handlebars::new();
    h.register_helper("format", Box::new(format));
    h.register_helper("fgt", Box::new(crate::fgt));

    let data = molecule_to_template_data(mol);
    h.render_template(template, &data).map_err(failure::err_msg)
}

#[test]
fn test_template_render() {
    let mols = io::read("tests/files/mol2/LTL-crysin-ds.mol2").expect("molecules");
    let mol = &mols[0];

    let template = io::read_file("tests/files/templates/xyz.hbs").expect("template xyz.hbs");
    let x = render_molecule_with(&mol, &template).unwrap();
}
// core:1 ends here
