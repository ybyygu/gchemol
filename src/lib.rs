// header
// #+name: bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d

// [[file:~/Workspace/Programming/gchemol/gchemol.note::bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d][bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d]]
//===============================================================================#
//   DESCRIPTION:  gchemol: a Graph-based CHEMical Objects Library
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  rewrite from scratch using rust
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 2 or upper
//       CREATED:  <2018-04-10 Tue 15:46>
//       UPDATED:  <2018-10-25 Thu 18:29>
//===============================================================================#
// bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d ends here

// base

// [[file:~/Workspace/Programming/gchemol/gchemol.note::*base][base:1]]
extern crate gchemol_geometry;
pub mod geometry {
    pub use gchemol_geometry::*;
}

pub use gchemol_core::*;

pub mod io {
    pub use gchemol_readwrite::io::*;
}

pub mod prelude {
    pub use gchemol_geometry::*;
    pub use gchemol_geometry::prelude::*;
    pub use gchemol_readwrite::prelude::*;
}
// base:1 ends here
