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
    /// Cartesian coordinates
    pub position: Point3D,
    /// Atom type, could be an element or a pseudo-atom
    pub kind: AtomKind,
    /// Atom nick name
    pub name: String,
    /// Would be managed by its parent molecule
    pub index: usize,
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            kind: Element(6),   // carbon atom
            position: [0.0; 3],
            name: "carbon".to_string(),
            index: 0,
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
// 150189fd-57d9-4e19-a888-d64497f5ba7e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::7e463bf4-a6ab-4648-a8a2-4b2023d1c588][7e463bf4-a6ab-4648-a8a2-4b2023d1c588]]
use std::str::FromStr;
use std::fmt;

use errors::*;

impl FromStr for Atom {
    type Err = Error;

    fn from_str(line: &str) -> Result<Self> {
        let parts: Vec<_> = line.split_whitespace().collect();
        if parts.len() != 4 {
            bail!("Incorrect number of data fields: {:?}", line);
        }

        let sym = parts[0];
        let msg = format!("Incorrect coordindate fields: {:}", parts[1]);
        let px: f64 = parts[1].parse().chain_err(|| msg)?;
        let msg = format!("Incorrect coordindate fields: {:}", parts[2]);
        let py: f64 = parts[2].parse().chain_err(|| msg)?;
        let msg = format!("Incorrect coordindate fields: {:}", parts[3]);
        let pz: f64 = parts[3].parse().chain_err(|| msg)?;

        let mut atom = Atom::new(sym, [px, py, pz]);
        atom.name = sym.to_string();

        Ok(atom)
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:6} {:-12.6} {:-12.6} {:-12.6}",
               self.symbol(),
               self.position[0],
               self.position[1],
               self.position[2]
        )
    }
}
// 7e463bf4-a6ab-4648-a8a2-4b2023d1c588 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b88435fd-d51c-48b8-880c-425b94b905e9][b88435fd-d51c-48b8-880c-425b94b905e9]]
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
}
// b88435fd-d51c-48b8-880c-425b94b905e9 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::cfdf0fc1-97a2-4da4-b0bb-a9baee31d275][cfdf0fc1-97a2-4da4-b0bb-a9baee31d275]]
#[test]
fn test_atom_string_conversion() {
    let line = "H 1.0 1.0 1.0";
    let a: Atom = line.parse().unwrap();
    assert_eq!(1, a.number());
    let line = a.to_string();
    let b: Atom = line.parse().unwrap();
    assert_eq!(a.kind, b.kind);
    assert_eq!(a.position, b.position);
    let line = "24 0.124 1.230 2.349";
    let a: Atom = line.parse().unwrap();
}
// cfdf0fc1-97a2-4da4-b0bb-a9baee31d275 ends here
