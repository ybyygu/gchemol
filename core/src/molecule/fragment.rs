// fragment.rs
// :PROPERTIES:
// :header-args: :tangle src/molecule/fragment.rs
// :END:
// Split molecule into fragments based on bonding connectivity.

// #+name: f0258648-03f4-41c9-949e-f3677c3b44bc

use super::*;

impl Molecule {
    /// Break molecule into smaller fragments based on its bonding connectivity
    pub fn fragment(&self) -> Vec<Molecule> {
        let graph = &self.graph;

        let mut mols = vec![];
        let subgraphs = connected_component_subgraphs(graph);

        for g in subgraphs {
            mols.push(Molecule::from_graph(g));
        }

        mols
    }

    /// Return a shallow connectivity graph without copying atom/bond data
    pub fn bond_graph(&self) -> UnGraph<usize, usize> {
        let graph = &self.graph;
        let mut result_g = UnGraph::with_capacity(graph.node_count(), graph.edge_count());
        // mapping from old node index to new node index
        let mut node_index_map = vec![AtomIndex::new(0); graph.node_count()];

        for (i, a) in self.view_atoms() {
            node_index_map[i-1] = result_g.add_node(i);
        }
        for (i, j, b) in self.view_bonds() {
            let source = node_index_map[i-1];
            let target = node_index_map[j-1];
            result_g.add_edge(source, target, 1);
        }
        result_g
    }

    /// Return the number of fragments in molecule
    pub fn nfragments(&self) -> usize {
        let gr = &self.bond_graph();
        // create a temporary graph
        let n = self.natoms();
        if n <= 1 {
            n
        } else {
            petgraph::algo::connected_components(&gr)
        }
    }
}


/// Generate connected components as subgraphs
fn connected_component_subgraphs(graph: &MolGraph) -> Vec<MolGraph>{
    // get fragments from connected components
    let components = petgraph::algo::kosaraju_scc(graph);

    let mut graphs = vec![];
    for nodes in components {
        let g: MolGraph = graph.filter_map
            (
                // node closure:
                // keep atoms in the same component
                |i, n| {
                    if nodes.contains(&i) {
                        Some(n.clone())
                    } else {
                        None
                    }
                },
                // edge closure:
                // keep the bond if bonded atoms are both in the same component
                |i, e| {
                    let (n1, n2) = graph.edge_endpoints(i).unwrap();
                    if nodes.contains(&n1) && nodes.contains(&n2) {
                        Some(e.clone())
                    } else {
                        None
                    }
                }
            );
        graphs.push(g);
    }

    graphs
}

// [[file:~/Workspace/Programming/gchemol/core/gchemol-core.note::9b65e04a-b7a8-4118-a2df-d0345c13832c][9b65e04a-b7a8-4118-a2df-d0345c13832c]]
impl Molecule {
    // FIXME: add bonds
    /// Create molecule from small fragments (molecules)
    pub fn combined(mols: &Vec<Molecule>) -> Self {
        let mut mol = Molecule::new("combined");
        for m in mols {
            // add atoms
            for a in m.atoms() {
                let n = mol.add_atom(a.clone());
            }
            // add bonds
            for b in m.bonds() {
                //
            }
        }

        mol
    }
}
// 9b65e04a-b7a8-4118-a2df-d0345c13832c ends here
