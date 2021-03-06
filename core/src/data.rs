// header

//===============================================================================#
//   DESCRIPTION:  provide important chemical data for atoms and molecules
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-12 Thu 14:40>
//       UPDATED:  <2019-07-03 Wed 16:59>
//===============================================================================#

use crate::molecule::{
    Atom,
    Bond,
    Molecule,
    AtomIndex,
    BondKind,
};

use crate::geometry::prelude::euclidean_distance;
use crate::core_utils::*;

// base

#[derive(Deserialize, Debug)]
#[serde(default)]
pub struct Config {
    bonding_ratio: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            bonding_ratio: 1.15,
        }
    }
}

// element radii data
// Element radii data taking from: https://mendeleev.readthedocs.io/en/stable/data.html

// Data in columns:
// covalent_radii_single covalent_radii_double, covalent_radii_triple, vdw_radii


const RADII_DATA: [[f64; 4]; 118] =
    [[0.32, 0.32, 0.32, 1.1],
     [0.46, 0.46, 0.46, 1.4],
     [1.33, 1.24, 1.24, 1.82],
     [1.02, 0.9, 0.85, 1.53],
     [0.85, 0.78, 0.73, 1.92],
     [0.75, 0.67, 0.6, 1.7],
     [0.71, 0.6, 0.54, 1.55],
     [0.63, 0.57, 0.53, 1.52],
     [0.64, 0.59, 0.53, 1.47],
     [0.67, 0.96, 0.96, 1.54],
     [1.55, 1.6, 1.6, 2.27],
     [1.39, 1.32, 1.27, 1.73],
     [1.26, 1.13, 1.11, 1.84],
     [1.16, 1.07, 1.02, 2.1],
     [1.11, 1.02, 0.94, 1.8],
     [1.03, 0.94, 0.95, 1.8],
     [0.99, 0.95, 0.93, 1.75],
     [0.96, 1.07, 0.96, 1.88],
     [1.96, 1.93, 1.93, 2.75],
     [1.71, 1.47, 1.33, 2.31],
     [1.48, 1.16, 1.14, 2.15],
     [1.36, 1.17, 1.08, 2.11],
     [1.34, 1.12, 1.06, 2.07],
     [1.22, 1.11, 1.03, 2.06],
     [1.19, 1.05, 1.03, 2.05],
     [1.16, 1.09, 1.02, 2.04],
     [1.11, 1.03, 0.96, 2.0],
     [1.1, 1.01, 1.01, 1.97],
     [1.12, 1.15, 1.2, 1.96],
     [1.18, 1.2, 1.2, 2.01],
     [1.24, 1.17, 1.21, 1.87],
     [1.21, 1.11, 1.14, 2.11],
     [1.21, 1.14, 1.06, 1.85],
     [1.16, 1.07, 1.07, 1.9],
     [1.14, 1.09, 1.1, 1.85],
     [1.17, 1.21, 1.08, 2.02],
     [2.1, 2.02, 2.02, 3.03],
     [1.85, 1.57, 1.39, 2.49],
     [1.63, 1.3, 1.24, 2.32],
     [1.54, 1.27, 1.21, 2.23],
     [1.47, 1.25, 1.16, 2.18],
     [1.38, 1.21, 1.13, 2.17],
     [1.28, 1.2, 1.1, 2.16],
     [1.25, 1.14, 1.03, 2.13],
     [1.25, 1.1, 1.06, 2.1],
     [1.2, 1.17, 1.12, 2.1],
     [1.28, 1.39, 1.37, 2.11],
     [1.36, 1.44, 1.44, 2.18],
     [1.42, 1.36, 1.46, 1.93],
     [1.4, 1.3, 1.32, 2.17],
     [1.4, 1.33, 1.27, 2.06],
     [1.36, 1.28, 1.21, 2.06],
     [1.33, 1.29, 1.25, 1.98],
     [1.31, 1.35, 1.22, 2.16],
     [2.32, 2.09, 2.09, 3.43],
     [1.96, 1.61, 1.49, 2.68],
     [1.8, 1.39, 1.39, 2.43],
     [1.63, 1.37, 1.31, 2.42],
     [1.76, 1.38, 1.28, 2.4],
     [1.74, 1.37, 1.37, 2.39],
     [1.73, 1.35, 1.35, 2.38],
     [1.72, 1.34, 1.34, 2.36],
     [1.68, 1.34, 1.34, 2.35],
     [1.69, 1.35, 1.32, 2.34],
     [1.68, 1.35, 1.35, 2.33],
     [1.67, 1.33, 1.33, 2.31],
     [1.66, 1.33, 1.33, 2.3],
     [1.65, 1.33, 1.33, 2.29],
     [1.64, 1.31, 1.31, 2.27],
     [1.7, 1.29, 1.29, 2.26],
     [1.62, 1.31, 1.31, 2.24],
     [1.52, 1.28, 1.22, 2.23],
     [1.46, 1.26, 1.19, 2.22],
     [1.37, 1.2, 1.15, 2.18],
     [1.31, 1.19, 1.1, 2.16],
     [1.29, 1.16, 1.09, 2.16],
     [1.22, 1.15, 1.07, 2.13],
     [1.23, 1.12, 1.1, 2.13],
     [1.24, 1.21, 1.23, 2.14],
     [1.33, 1.42, 1.42, 2.23],
     [1.44, 1.42, 1.5, 1.96],
     [1.44, 1.35, 1.37, 2.02],
     [1.51, 1.41, 1.35, 2.07],
     [1.45, 1.35, 1.29, 1.97],
     [1.47, 1.38, 1.38, 2.02],
     [1.42, 1.45, 1.33, 2.2],
     [2.23, 2.18, 2.18, 3.48],
     [2.01, 1.73, 1.59, 2.83],
     [1.86, 1.53, 1.4, 2.47],
     [1.75, 1.43, 1.36, 2.45],
     [1.69, 1.38, 1.29, 2.43],
     [1.7, 1.34, 1.18, 2.41],
     [1.71, 1.36, 1.16, 2.39],
     [1.72, 1.35, 1.35, 2.43],
     [1.66, 1.35, 1.35, 2.44],
     [1.66, 1.36, 1.36, 2.45],
     [1.68, 1.39, 1.39, 2.44],
     [1.68, 1.4, 1.4, 2.45],
     [1.65, 1.4, 1.4, 2.45],
     [1.67, 1.67, 1.67, 2.45],
     [1.73, 1.39, 1.39, 2.46],
     [1.76, 1.76, 1.76, 2.46],
     [1.61, 1.41, 1.41, 2.46],
     [1.57, 1.4, 1.31, 2.46],
     [1.49, 1.36, 1.26, 2.46],
     [1.43, 1.28, 1.21, 2.46],
     [1.41, 1.28, 1.19, 2.46],
     [1.34, 1.25, 1.18, 2.46],
     [1.29, 1.25, 1.13, 2.46],
     [1.28, 1.16, 1.18, 2.46],
     [1.21, 1.16, 1.18, 2.46],
     [1.22, 1.37, 1.3, 2.46],
     [1.36, 1.36, 1.36, 2.46],
     [1.43, 1.43, 1.43, 2.46],
     [1.62, 1.62, 1.62, 2.46],
     [1.75, 1.75, 1.75, 2.46],
     [1.65, 1.65, 1.65, 2.46],
     [1.57, 1.57, 1.57, 2.46]
    ];

// methods for Atom

/// Return covalent radius for single, double, or triple bonds
fn get_cov_radius(element_number: usize, bond_order: usize) -> Option<f64> {
    if element_number <= RADII_DATA.len() {
        // only for single, double, or triple bond
        if bond_order > 0 && bond_order <= 3 {
            let irow = element_number - 1;
            let icol = bond_order - 1;

            let r = RADII_DATA[irow][icol];
            return Some(r);
        }
    }

    None
}

/// Return Van der Waals radius if any
fn get_vdw_radius(element_number: usize) -> Option<f64> {
    // column index to vdw radii
    let icol = 3;
    if element_number <= RADII_DATA.len() {
        let irow = element_number - 1;
        let r = RADII_DATA[irow][icol];
        return Some(r);
    }
    None
}

pub fn guess_bond_kind(atom1: &Atom, atom2: &Atom) -> BondKind {
    let config = match envy::prefixed("GCHEMOL_").from_env::<Config>() {
        Ok(config) => config,
        Err(error) => {
            error!("{:?}", error);
            Config::default()
        }
    };

    if let Some(cr1) = get_cov_radius(atom1.number(), 1) {
        if let Some(cr2) = get_cov_radius(atom2.number(), 1) {
            let d12 = atom1.distance(atom2);

            let r = config.bonding_ratio;
            let rcutoff = (cr1 + cr2) * r;
            if d12 > rcutoff {
                return BondKind::Dummy;
            } else {
                return BondKind::Single;
            }
        }
    }

    BondKind::Dummy
}

impl Atom {
    /// Access covalent atomic radii for single bond
    /// Return None if no data available
    pub fn cov_radius(&self) -> Option<f64> {
        get_cov_radius(self.number(), 1)
    }

    pub fn vdw_radius(&self) -> Option<f64> {
        get_vdw_radius(self.number())
    }
}

// src
// to be removed

impl Molecule {
    /// Deduce the distance bound for given pair of atoms i and j according
    /// their bonding pattern
    pub fn distance_bound(&self, i: AtomIndex, j: AtomIndex) -> Result<[f64; 2]> {
        let min_rij = 0.5;
        let max_rij = 40.0;

        let atom_i = self.get_atom(i).ok_or(format_err!("no atom {:?} in molecule", i))?;
        let atom_j = self.get_atom(j).ok_or(format_err!("no atom {:?} in molecule", j))?;

        // 1. get data
        // vdw radii as cutoff for non-bonding interaction
        let rwi = atom_i.vdw_radius().ok_or(format_err!("no radius data"))?;
        let rwj = atom_j.vdw_radius().ok_or(format_err!("no radius data"))?;
        let srwij = rwi + rwj;

        // covalent radii as cutoff for bonding interaction
        let rci1 = atom_i.cov_radius().ok_or(format_err!("no radius data"))?;
        let rcj1 = atom_j.cov_radius().ok_or(format_err!("no radius data"))?;
        // radii for double bond
        let rci2 = get_cov_radius(atom_i.number(), 2).ok_or(format_err!("no radius data"))?;
        let rcj2 = get_cov_radius(atom_j.number(), 2).ok_or(format_err!("no radius data"))?;
        // radii for triple bond
        let rci3 = get_cov_radius(atom_i.number(), 3).ok_or(format_err!("no radius data"))?;
        let rcj3 = get_cov_radius(atom_j.number(), 3).ok_or(format_err!("no radius data"))?;

        let srcij1 = rci1 + rcj1;
        let srcij2 = rci2 + rcj2;
        let srcij3 = rci3 + rcj3;

        let mut bound = [min_rij, max_rij];
        // 2. check bonding
        if let Some(bond) = self.get_bond_between(i, j) {
            match bond.kind {
                // treat dummy bond as non-bonding
                BondKind::Dummy => {
                    bound = [srwij, max_rij];
                },
                BondKind::Partial => {
                    bound = [srcij1 * 1.2, srwij];
                },
                BondKind::Single => {
                    let lij = srcij1 * 0.8;
                    let uij = srcij1 * 1.2;
                    bound = [lij, uij];
                },
                BondKind::Aromatic => {
                    bound = [srcij2, srcij1];
                },
                BondKind::Double => {
                    bound = [srcij2 * 0.8, srcij1];
                },
                BondKind::Triple => {
                    bound = [srcij3 * 0.8, srcij2];
                },
                BondKind::Quadruple => {
                    bound = [srcij3*0.6, srcij3];
                },
            }
        } else {
            // non-bonding
            // treat as weak bonding if atom pair are already close to each other
            let dij = atom_i.distance(atom_j);
            if dij > srcij1*1.2 && dij < srwij {
                bound = [srcij1*1.2, srwij];
            } else {
                bound = [srwij, max_rij];
            }
        }

        // 3. final check
        if bound[0] > bound[1] {
            bound.swap(0, 1);
        }

        Ok(bound)
    }
}

// test

#[test]
fn test_guess_bond_kind() {
    // CH4 molecule
    let atom1 = Atom::new("C", [-0.90203687,  0.62555259,  0.0081889 ]);
    let atom2 = Atom::new("H", [-0.54538244, -0.38325741,  0.0081889 ]);
    let atom3 = Atom::new("H", [-0.54536403,  1.12995078,  0.88184041]);
    let atom4 = Atom::new("H", [-0.54536403,  1.12995078, -0.8654626 ]);
    let atom5 = Atom::new("H", [-1.97203687,  0.62556577,  0.0081889 ]);

    let x = guess_bond_kind(&atom1, &atom2);
    assert!(x != BondKind::Dummy);
    let x = guess_bond_kind(&atom2, &atom3);
    assert!(x == BondKind::Dummy);

    assert!(atom1.cov_radius().is_some());
    assert!(atom1.vdw_radius().is_some());
}
