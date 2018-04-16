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
    pub kind : BondKind,
    pub name : String,
    /// set this attribute for arbitrary bond order
    order    : Option<f64>,
}

impl Default for Bond {
    fn default() -> Self {
        Bond {
            order : None,
            kind  : BondKind::Single,
            name  : String::default(),
        }
    }
}

impl Bond {
    /// Return bond order
    pub fn order(&self) -> f64 {
        if let Some(order) = self.order {
            order
        } else {
            match self.kind {
                BondKind::Dummy     => 0.0,
                BondKind::Partial   => 0.5,
                BondKind::Single    => 1.0,
                BondKind::Aromatic  => 1.5,
                BondKind::Double    => 2.0,
                BondKind::Triple    => 3.0,
                BondKind::Quadruple => 4.0,
            }
        }
    }
}
// 75287fd8-649d-496d-8c50-40f9247a4c10 ends here
