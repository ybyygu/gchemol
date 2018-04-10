// [[file:~/Workspace/Programming/gchemol/gchemol.note::942dedaa-9351-426e-9be9-cdb640ec2b75][942dedaa-9351-426e-9be9-cdb640ec2b75]]
use petgraph::prelude::*;
use petgraph::graph::edge_index;

use Point3D;
use Points;
use Atom;

type AtomIndex = NodeIndex<usize>;
type BondIndex = EdgeIndex<usize>;

pub struct Molecule {
    pub atoms: Vec<Atom>,
    graph: StableUnGraph<Atom, usize>,
}

impl Molecule {
    /// get positions of all atoms
    fn positions(&self) -> impl Iterator<Item = &Point3D> {
        self.atoms.iter().map(|ref a| &a.position)
    }

    fn symbols(&self) -> impl Iterator<Item = &str> {
        self.atoms.iter().map(|a| a.symbol.as_str())
    }
}

#[test]
fn test_graph() {
    let mut g = StableUnGraph::<_, _>::default();
    let a = g.add_node(Atom::default());
    let b = g.add_node(Atom::default());
    let c = g.add_node(Atom::default());
    let d = g.add_node(Atom::default());
    g.extend_with_edges(&[
        (a, b, 1),
        (a, c, 2),
        (b, c, 3),
        (c, c, 4),
        (a, d, 5),
    ]);

    println!("{:?}", g[a]);
    println!("{:?}", g[a]);
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here
