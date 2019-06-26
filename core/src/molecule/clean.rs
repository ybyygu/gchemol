// bounds

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*bounds][bounds:1]]
use super::*;

type Bounds = HashMap<(AtomIndex, AtomIndex), f64>;

// return distance bounds between atoms
// upper-tri for upper bounds
// lower-tri for lower bounds
fn get_distance_bounds_v1(mol: &Molecule) -> Bounds {
    // max distance between two atoms
    let max_rij = 90.0;

    let mut bounds = HashMap::new();
    let node_indices = mol.sites();
    let nnodes = node_indices.len();

    for i in 0..nnodes {
        for j in (i + 1)..nnodes {
            let node_i = node_indices[i];
            let node_j = node_indices[j];
            let atom_i = &mol.graph[node_i];
            let atom_j = &mol.graph[node_j];

            // use vdw radii as the lower bound for non-bonded pair
            let vri = atom_i.vdw_radius().unwrap();
            let vrj = atom_j.vdw_radius().unwrap();
            let vrij = vri + vrj;

            // use covalent radii as the lower bound for bonded pair
            let cri = atom_i.cov_radius().unwrap();
            let crj = atom_j.cov_radius().unwrap();
            let crij = cri + crj;

            let lij = crij * 0.8;
            let uij = vrij * 0.8;
            let uij = if uij > crij * 1.2 { uij } else { crij * 1.2 };

            debug_assert!(lij <= uij);
            let dij = atom_i.distance(atom_j);
            // if i and j is directly bonded
            // set covalent radius as the lower bound
            // or set vdw radius as the lower bound if not bonded
            if let Some(nb) = mol.nbonds_between(node_i, node_j) {
                if nb == 1 {
                    if dij >= lij && dij < crij * 1.2 {
                        bounds.insert((node_i, node_j), dij);
                        bounds.insert((node_j, node_i), dij);
                    } else {
                        bounds.insert((node_i, node_j), lij);
                        bounds.insert((node_j, node_i), lij);
                    }
                } else if nb == 2 {
                    if dij > uij && dij < max_rij {
                        bounds.insert((node_i, node_j), dij);
                        bounds.insert((node_j, node_i), dij);
                    } else {
                        bounds.insert((node_i, node_j), uij);
                        bounds.insert((node_j, node_i), uij + dij);
                    }
                } else {
                    if dij > uij && dij < max_rij {
                        bounds.insert((node_i, node_j), dij);
                        bounds.insert((node_j, node_i), max_rij);
                    } else {
                        bounds.insert((node_i, node_j), uij);
                        bounds.insert((node_j, node_i), max_rij);
                    }
                }
            } else {
                bounds.insert((node_i, node_j), uij);
                bounds.insert((node_j, node_i), max_rij);
            }
            // println!("pair: {:3}-{:3}, lij = {:-4.1}, uij = {:-4.1}", node_i.index() + 1, node_j.index() + 1, lij, uij);
        }
    }

    bounds
}
// bounds:1 ends here

// core

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*core][core:1]]
// the weight between two atoms
fn get_weight_between(lij: f64, uij: f64, dij: f64) -> f64 {
    debug_assert!(lij <= uij);

    let weight = if dij >= lij && dij < uij {
        // avoid dividing by zero (nan)
        1e-4
    } else if dij < lij {
        1.0
    } else {
        1.0
    };

    weight
}

impl Molecule {
    pub fn set_momentum(&mut self, index: AtomIndex, m: Point3D) {
        let mut atom = &mut self.graph[index];
        atom.set_momentum(m);
    }

    /// Clean up molecule geometry using stress majorization algorithm
    pub fn clean(&mut self) -> Result<()> {
        let bounds = get_distance_bounds_v1(&self);
        // fix_13_distance(&mut bounds, &self);
        // let bounds = get_bounds();
        // print_bounds(&bounds, self.natoms());
        let node_indices: Vec<_> = self.graph.node_indices().collect();
        let nnodes = node_indices.len();

        let maxcycle = nnodes * 100;
        let mut icycle = 0;
        let ecut = 1E-4;

        let mut old_stress = 0.0;
        loop {
            let mut stress = 0.0;
            let mut positions_new = vec![];
            for i in 0..nnodes {
                let node_i = node_indices[i];
                let mut pi = self
                    .get_atom(node_i)
                    .expect("atom i from node_i")
                    .position();
                let mut pi_new = [0.0; 3];
                let mut wijs = vec![];
                let npairs = (nnodes - 1) as f64;
                let mut stress_i = 0.0;
                for j in 0..nnodes {
                    // skip self-interaction
                    if i == j {
                        continue;
                    };

                    let node_j = node_indices[j];
                    let pj = self
                        .get_atom(node_j)
                        .expect("atom j from node_j")
                        .position();

                    // current distance
                    let cur_dij = euclidean_distance(pi, pj);

                    // lower bound and upper bound for pair distance
                    let lij = if node_i < node_j {
                        bounds[&(node_i, node_j)]
                    } else {
                        bounds[&(node_j, node_i)]
                    };
                    // let uij = bounds[&(node_j, node_i)];
                    let uij = if node_i < node_j {
                        bounds[&(node_j, node_i)]
                    } else {
                        bounds[&(node_i, node_j)]
                    };
                    // let (lij, uij) = (lij.min(uij), uij.max(lij));

                    // ij pair counts twice, so divide the weight
                    // let wij = lij.powi(-4);
                    let wij = get_weight_between(lij, uij, cur_dij);
                    // dbg!((wij, lij, uij, cur_dij));
                    let wij = 0.5 * wij;
                    wijs.push(wij);

                    // collect position contribution of atom j to atom i
                    let xij = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
                    let mut fij = [0.0; 3];
                    for v in 0..3 {
                        pi_new[v] += wij * (pj[v] + lij / cur_dij * (pi[v] - pj[v]));
                        fij[v] = (1.0 - lij / cur_dij) * xij[v];
                    }
                    stress_i += wij * (cur_dij - lij).powi(2);
                }

                // weight sum
                let swij: f64 = wijs.iter().sum();

                // FIXME: if all pair weights are zero
                // dbg!(swij);
                debug_assert!(swij.abs() >= 1e-4);
                for v in 0..3 {
                    pi_new[v] /= swij;
                }
                positions_new.push((node_i, pi_new));
                stress += stress_i;
            }

            debug!("cycle: {} energy = {:?}", icycle, stress);
            // println!("cycle: {} energy = {:?}", icycle, stress);

            // update positions
            for (node, position) in positions_new {
                self.set_position(node, position);
            }

            if stress.is_nan() {
                bail!("found invalid number: {:?}", stress);
            }

            if stress < ecut
                || (stress - old_stress).abs() < ecut
                || (stress - old_stress).abs() / stress < ecut
            {
                break;
            }

            icycle += 1;
            if icycle > maxcycle {
                break;
            }

            old_stress = stress;
        }

        Ok(())
    }
}
// core:1 ends here
