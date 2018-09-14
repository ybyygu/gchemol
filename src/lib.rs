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
//       UPDATED:  <2018-09-14 Fri 18:59>
//===============================================================================#
// bdab2ff7-59d6-4b5e-8b47-53eaccf5e64d ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::53cbd3c0-e164-4bad-b535-6fd6df916650][53cbd3c0-e164-4bad-b535-6fd6df916650]]
extern crate gchemol_geometry;
pub mod geometry {
    pub use gchemol_geometry::*;
}

extern crate gchemol_core;
pub use gchemol_core::*;

extern crate gchemol_readwrite;
pub mod io {
    pub use gchemol_readwrite::io::*;
}

pub mod prelude {
    pub use gchemol_readwrite::prelude::*;
    pub use gchemol_geometry::*;
}
// 53cbd3c0-e164-4bad-b535-6fd6df916650 ends here
