// [[file:~/Workspace/Programming/gchemol/readwrite/readwrite.note::abd992a3-0aff-40a9-ab09-7d4956ce57ee][abd992a3-0aff-40a9-ab09-7d4956ce57ee]]
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
mod parser;
// abd992a3-0aff-40a9-ab09-7d4956ce57ee ends here
