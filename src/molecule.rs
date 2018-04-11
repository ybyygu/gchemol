// [[file:~/Workspace/Programming/gchemol/gchemol.note::942dedaa-9351-426e-9be9-cdb640ec2b75][942dedaa-9351-426e-9be9-cdb640ec2b75]]
use Point3D;
use Points;
use Atom;

pub struct Molecule {
    pub atoms: Vec<Atom>,
}

impl Molecule {
    fn new() -> Self {
        Molecule {
            atoms: vec![],
        }
    }

    /// get positions of all atoms
    fn positions(&self) -> impl Iterator<Item = &Point3D> {
        self.atoms.iter().map(|ref a| &a.position)
    }

    fn symbols(&self) -> impl Iterator<Item = &str> {
        self.atoms.iter().map(|ref a| a.symbol())
    }

    /// Add a single atom into molecule
    fn add_atom(&mut self, atom: Atom) {
        // let index = self.add_node(atom);
    }
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here
