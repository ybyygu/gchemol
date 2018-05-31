// [[file:~/Workspace/Programming/gchemol/gchemol.note::6b59958f-7e56-4b16-b02a-cc01e5de3da8][6b59958f-7e56-4b16-b02a-cc01e5de3da8]]
#[macro_use]
extern crate serde_json;
#[macro_use]
extern crate serde_derive;

extern crate gchemol;
extern crate handlebars;


use gchemol::errors::*;
use handlebars::{
    to_json,
    Handlebars,
    Helper,
    HelperResult,
    RenderContext,
    RenderError
};

static TEMPLATE: &'static str =
"{{natoms}}
{{title}}
{{#each atoms as |a| ~}}
{{format a.symbol width=3 align=\"center\"}} {{format a.x width=6 prec=2}} {{format a.y}} {{format a.z}}
{{/each~}}
";
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
            .ok_or(RenderError::new("Wrong type of param 0"))?;
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
            .ok_or(RenderError::new("Wrong type of param 0"))?;

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
        // TODO
    }

    Ok(())
}

// // define a custom helper for formatting float number
// fn format_number(h: &Helper, _: &Handlebars, rc: &mut RenderContext) -> std::result::Result<(), RenderError> {
//     // get parameter from helper or throw an error
//     let param = h.param(0).ok_or(RenderError::new("Param 0 is required for format helper."))?;
//     let num: f64 = param.value().as_f64().ok_or(RenderError::new("Param 0 sould be a float number"))?;
//     let rendered = format!("{:-8.4}", num);
//     rc.writer.write(rendered.into_bytes().as_ref())?;
//     Ok(())
// }

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
    let mols = gchemol::io::read("tests/files/mol2/LTL-crysin-ds.mol2")?;
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

    let mut h = Handlebars::new();
    h.register_helper("format", Box::new(format));
    let x = h.render_template(TEMPLATE, &data).chain_err(|| "failed")?;
    println!("{:}", x);

    Ok(())
}
// 6b59958f-7e56-4b16-b02a-cc01e5de3da8 ends here
