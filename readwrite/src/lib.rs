// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

#[macro_use] extern crate nom;
#[macro_use] extern crate quicli;
#[macro_use] extern crate handlebars;
#[macro_use] extern crate indexmap;
#[macro_use] extern crate serde_json;
#[macro_use] extern crate serde_derive;

// for tests only
#[cfg(test)]
#[macro_use] extern crate approx;

extern crate serde;
extern crate itertools;

extern crate gchemol_core;

#[macro_use] pub mod parser;

pub mod io;
pub mod formats;
pub mod template;

pub mod prelude {
    pub use io::prelude::*;
    pub use template::*;
}
