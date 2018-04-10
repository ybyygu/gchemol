// [[file:~/Workspace/Programming/gchemol/gchemol.note::150189fd-57d9-4e19-a888-d64497f5ba7e][150189fd-57d9-4e19-a888-d64497f5ba7e]]
use std::hash::{Hash, Hasher};
use std::cmp::Ordering;

use Point3D;

#[derive (Debug, Clone)]
/// simple atom data structure
pub struct Atom {
    pub symbol: String,
    pub position: Point3D,
}

impl Atom {
    pub fn new<T: Into<String>>(symbol: T, position: Point3D) -> Self {
        Atom {
            symbol: symbol.into(),
            position: position,
        }
    }
}

impl Default for Atom {
    fn default() -> Self {
        Atom::new("C", [0.0; 3])
    }
}
// 150189fd-57d9-4e19-a888-d64497f5ba7e ends here
