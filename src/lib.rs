// header
// #+name: bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d

//===============================================================================#
//   DESCRIPTION:  gchemol: a Graph-based CHEMical Objects Library
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  rewrite from scratch using rust
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 2 or upper
//       CREATED:  <2018-04-10 Tue 15:46>
//       UPDATED:  <2018-09-21 Fri 10:10>
//===============================================================================#

// base
// #+name: 53cbd3c0-e164-4bad-b535-6fd6df916650

extern crate gchemol_readwrite;
#[macro_use]  pub mod parser {
    pub use gchemol_readwrite::parser::*;
}

extern crate gchemol_geometry;
pub mod geometry {
    pub use gchemol_geometry::*;
}

extern crate gchemol_core;
pub use gchemol_core::*;

pub mod io {
    pub use gchemol_readwrite::io::*;
}

pub mod prelude {
    pub use gchemol_readwrite::prelude::*;
    pub use gchemol_geometry::*;
}
