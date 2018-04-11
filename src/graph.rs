// [[file:~/Workspace/Programming/gchemol/gchemol.note::cffd4016-440e-47bb-8660-0a421e930e32][cffd4016-440e-47bb-8660-0a421e930e32]]
use std::collections::HashMap;

use petgraph::prelude::*;
use petgraph as pg;

use Point3D;
use Points;
use Atom;
// cffd4016-440e-47bb-8660-0a421e930e32 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::94a2f00a-b902-4997-9ba7-02f354fbaee2][94a2f00a-b902-4997-9ba7-02f354fbaee2]]
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
    println!("{:?}", a);
}
// 94a2f00a-b902-4997-9ba7-02f354fbaee2 ends here
