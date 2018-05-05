// [[file:~/Workspace/Programming/gchemol/gchemol.note::d333cb1f-e622-462f-a892-4906c85b7da0][d333cb1f-e622-462f-a892-4906c85b7da0]]
#[derive(Debug, Clone)]
pub struct AtomData {
    /// Atom type, could be an element or a pseudo-atom
    kind: AtomKind,
    /// Cartesian coordinates
    position: Point3D,
    /// Atom nick name
    name: String,
    /// Atomic momentum vector
    momentum: Point3D,
}

impl Default for AtomData {
    fn default()  -> Self {
        AtomData {
            kind: Element(6),   // carbon atom
            position: [0.0; 3],
            momentum: [0.0; 3],
            name: "carbon".into(),
        }
    }
}

impl AtomData {
    pub fn new() -> Self {
        AtomData::default()
    }

    /// Set atom position
    #[inline]
    pub fn position(&mut self, x: f64, y: f64, z: f64) -> &mut Self {
        self.position = [x, y, z]; self
    }

    /// Set atom kink using element number
    #[inline]
    pub fn element(&mut self, n: usize) -> &mut Self {
        self.kind = Element(n); self
    }

    /// Set atom kind using element symbol
    #[inline]
    pub fn symbol<T: Into<String>>(&mut self, s: T) -> &mut Self {
        self.kind = atom_kind_from_string(s.into()); self
    }

    #[inline]
    pub fn momentum(&mut self, x: f64, y: f64, z: f64) -> &mut Self {
        self.momentum = [x, y, z]; self
    }

    /// return a new `Atom` struct
    pub fn finish(&self) -> Atom {
        let mut atom = Atom::default();
        atom.data = self.clone();

        atom
    }
}

impl Atom {
    pub fn new<T: Into<String>>(s: T, p: Point3D) -> Self {
        AtomData::new()
            .symbol(s)
            .position(p[0], p[1], p[2])
            .finish()
    }

    pub fn build() -> AtomData {
        AtomData::new()
    }
}

#[test]
fn test_atom_builder() {
    let a = Atom::build()
        .position(0.0, 0.0, 1.2)
        .symbol("Fe")
        .element(13)
        .momentum(0.2, 0.2, 0.3)
        .finish();
}
// d333cb1f-e622-462f-a892-4906c85b7da0 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::150189fd-57d9-4e19-a888-d64497f5ba7e][150189fd-57d9-4e19-a888-d64497f5ba7e]]
use std::hash::{Hash, Hasher};
use std::cmp::Ordering;

use element::atom_kind_from_string;
use Point3D;
use AtomKind;
use Element;
use petgraph;

type AtomIndex = petgraph::prelude::NodeIndex;

#[derive (Debug, Clone)]
/// simple atom data structure
pub struct Atom {
    /// Would be managed by its parent molecule
    pub index: AtomIndex,

    /// private atom data
    data: AtomData,
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            index: AtomIndex::new(0),
            data: AtomData::default(),
        }
    }
}

impl Atom {
    /// shortcut for accessing atom symbol
    pub fn symbol(&self) -> &str {
        self.data.kind.symbol()
    }

    /// shortcut for accessing atom number
    pub fn number(&self) -> usize {
        self.data.kind.number()
    }

    /// shortcut for accessing atom position
    pub fn position(&self) -> Point3D {
        self.data.position
    }

    #[inline]
    pub fn set_position(&mut self, p: Point3D) {
        self.data.position = p;
    }

    #[inline]
    pub fn set_momentum(&mut self, m: Point3D) {
        self.data.momentum = m;
    }

    pub fn name(&self) -> &str {
        &self.data.name
    }

    #[inline]
    pub fn set_name<T: Into<String>>(&mut self, s: T) {
        self.data.name = s.into();
    }
}
// 150189fd-57d9-4e19-a888-d64497f5ba7e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b6d1e417-27da-4384-879a-db28960ed161][b6d1e417-27da-4384-879a-db28960ed161]]
use geometry::euclidean_distance;

impl Atom {
    pub fn distance(&self, other: &Atom) -> f64 {
        euclidean_distance(self.position(), other.position())
    }
}
// b6d1e417-27da-4384-879a-db28960ed161 ends here

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
        atom.set_name(sym);

        Ok(atom)
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:6} {:-12.6} {:-12.6} {:-12.6}",
               self.symbol(),
               self.data.position[0],
               self.data.position[1],
               self.data.position[2]
        )
    }
}
// 7e463bf4-a6ab-4648-a8a2-4b2023d1c588 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b88435fd-d51c-48b8-880c-425b94b905e9][b88435fd-d51c-48b8-880c-425b94b905e9]]
#[test]
fn test_atom_init() {
    let atom = Atom::default();
    let atom = Atom::new("Fe", [9.3; 3]);
    assert_eq!(9.3, atom.position()[0]);
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
    assert_eq!(a.symbol(), b.symbol());
    assert_eq!(a.position(), b.position());
    let line = "24 0.124 1.230 2.349";
    let a: Atom = line.parse().unwrap();
}
// cfdf0fc1-97a2-4da4-b0bb-a9baee31d275 ends here
