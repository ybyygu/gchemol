// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/readwrite/readwrite.note::*lib.rs][lib.rs:1]]
#[macro_use]
extern crate nom;
#[macro_use]
extern crate indexmap;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate handlebars;
#[macro_use]
extern crate serde_json;

// for tests only
#[cfg(test)]
#[macro_use]
extern crate approx;

extern crate gchemol_core;
extern crate itertools;
extern crate serde;

pub mod formats;
pub mod io;
pub mod template;

mod core_utils {
    pub use guts::prelude::*;
}

// FIXME: why must be located in lib.rs?
handlebars_helper!(fgt: |x: f64, y: f64| x > y);

pub mod prelude {
    pub use crate::io::prelude::*;
    pub use crate::template::*;
}
// lib.rs:1 ends here
