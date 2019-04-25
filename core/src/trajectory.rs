// imports

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*imports][imports:1]]
use crate::lattice::Lattice;
use crate::molecule::Molecule;

use crate::core_utils::*;
//use quicli::prelude::*;
// imports:1 ends here

// core

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*core][core:1]]
pub struct Conformation {}

pub struct Configuration {
    /// The molecule positions.
    positions: Vec<[f64; 3]>,

    lattice: Option<Lattice>,
}

pub struct Trajectory {
    /// Parent molecule of configuration frames in the trajectory.
    parent: Molecule,

    /// Frames represent molecular configuations along time axis.
    frames: Vec<Configuration>,
}

impl Trajectory {
    /// Construct `Trajectory` object from a list of `Molecule`
    pub fn new(mols: &[Molecule]) -> Self {
        let frames: Vec<_> = mols
            .iter()
            .map(|mol| Configuration {
                positions: mol.positions(),
                lattice: mol.lattice.clone(),
            })
            .collect();

        Self {
            frames,
            parent: mols[0].clone(),
        }
    }

    /// Return the number of frames in the trajectory.
    pub fn nframes(&self) -> usize {
        self.frames.len()
    }

    /// Return the number of atoms in each frame.
    pub fn natoms(&self) -> usize {
        self.parent.natoms()
    }
}
// core:1 ends here

// pub

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*pub][pub:1]]
use std::convert::TryFrom;

impl TryFrom<Vec<Molecule>> for Trajectory {
    type Error = Error;

    fn try_from(mols: Vec<Molecule>) -> Result<Self> {
        for (i, pair) in mols.windows(2).enumerate() {
            if !matchable(&pair[0], &pair[1]) {
                bail!("found inconsistent molecules: {} -- {}!", i, i + 1)
            }
        }

        // FIXME: avoid re-allocation
        Ok(Self::new(&mols))
    }
}

impl Trajectory {
    pub fn iter(&self) -> impl Iterator<Item = Molecule> {
        let mut mol = self.parent.clone();
        std::iter::from_fn(move || Some(mol.clone()))
    }
}
// pub:1 ends here

// conversion

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*conversion][conversion:1]]
// Check if two molecules are in different configuration.
fn matchable(mol1: &Molecule, mol2: &Molecule) -> bool {
    // check number of atoms
    let n1 = mol1.natoms();
    let n2 = mol2.natoms();
    if n1 != n2 {
        error!("different molecule size: {} != {}", n1, n2);
        return false;
    }

    // check elements
    let syms1 = mol1.symbols();
    let syms2 = mol2.symbols();

    for i in 0..n1 {
        let e1 = syms1[i];
        let e2 = syms2[i];
        if e1 != e2 {
            error!(
                "different element symbol for atom {}: {} != {}",
                i + 1,
                e1,
                e2
            );
            return false;
        }
    }

    // check lattice
    if mol1.lattice.is_none() && mol2.lattice.is_none() {
        true
    } else if mol1.lattice.is_some() && mol2.lattice.is_some() {
        true
    } else {
        error!("different lattice!");
        false
    }
}
// conversion:1 ends here

// interplate

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*interplate][interplate:1]]
impl Molecule {
    /// Interpolate between this structure and end_structure. Useful for
    /// construction of NEB inputs.
    fn interpolate(&self, end_mol: &Molecule, nimages: usize) -> Trajectory {
        unimplemented!()
    }
}
// interplate:1 ends here
