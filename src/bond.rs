// [[file:~/Workspace/Programming/gchemol/gchemol.note::75287fd8-649d-496d-8c50-40f9247a4c10][75287fd8-649d-496d-8c50-40f9247a4c10]]
use Atom;

/// https://en.wikipedia.org/wiki/Bond_order
#[derive(Debug, Clone, Copy, PartialOrd, PartialEq)]
pub enum BondKind {
    Dummy,
    Partial,
    Single,
    Aromatic,
    Double,
    Triple,
    Quadruple,
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub kind: BondKind,
    pub name: String,
}

impl Default for Bond {
    fn default() -> Self {
        Bond {
            kind: BondKind::Single,
            name: String::default(),
        }
    }
}
// 75287fd8-649d-496d-8c50-40f9247a4c10 ends here
