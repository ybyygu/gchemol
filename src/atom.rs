// [[file:~/Workspace/Programming/gchemol/gchemol.note::150189fd-57d9-4e19-a888-d64497f5ba7e][150189fd-57d9-4e19-a888-d64497f5ba7e]]
use std::hash::{Hash, Hasher};
use std::cmp::Ordering;

use element::atom_kind_from_string;
use Point3D;
use AtomKind;
use Element;

#[derive (Debug, Clone)]
/// simple atom data structure
pub struct Atom {
    /// Atom type, could be an element or a pseudo-atom
    pub kind: AtomKind,
    /// Cartesian coordinates
    pub position: Point3D,
    /// Atom nick name
    pub name: String,
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            kind: Element(6),   // carbon atom
            position: [0.0; 3],
            name: "carbon".to_string(),
        }
    }
}

impl Atom {
    pub fn new<T: Into<String>>(symbol: T, position: Point3D) -> Self {
        let element = atom_kind_from_string(symbol);
        let name = element.symbol().to_string();
        Atom {
            kind: element,
            name: name,
            position: position,
            ..Default::default()
        }
    }

    /// shortcut for accessing atom symbol
    pub fn symbol(&self) -> &str {
        self.kind.symbol()
    }

    /// shortcut for accessing atom number
    pub fn number(&self) -> usize {
        self.kind.number()
    }
}

#[test]
fn test_atom_init() {
    let atom = Atom::default();
    let atom = Atom::new("Fe", [9.3; 3]);
    assert_eq!(9.3, atom.position[0]);
    assert_eq!("Fe", atom.symbol());
    assert_eq!(26, atom.number());

    let atom = Atom::new("dummy", [9.3; 3]);
    assert_eq!("dummy", atom.symbol());
    assert_eq!(0, atom.number());
    println!("{:?}", atom.name);
}
// 150189fd-57d9-4e19-a888-d64497f5ba7e ends here
