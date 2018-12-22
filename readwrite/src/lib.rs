// lib.rs
// :PROPERTIES:
// :header-args: :tangle src/lib.rs
// :END:

// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::*lib.rs][lib.rs:1]]
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

pub mod io;
pub mod formats;
pub mod template;

mod core_utils {
    pub use quicli::prelude::*;
    pub type Result<T> = ::std::result::Result<T, Error>;
}

pub mod prelude {
    pub use crate::io::prelude::*;
    pub use crate::template::*;
}
// lib.rs:1 ends here
