// [[file:~/Workspace/Programming/gchemol/gchemol.note::75287fd8-649d-496d-8c50-40f9247a4c10][75287fd8-649d-496d-8c50-40f9247a4c10]]
use Atom;

#[derive(Debug, Clone)]
pub struct Bond {
    order: f64,
    name: String,
}

impl Default for Bond {
    fn default() -> Self {
        Bond {
            order: 1.0,
            name: String::default(),
        }
    }
}
// 75287fd8-649d-496d-8c50-40f9247a4c10 ends here
