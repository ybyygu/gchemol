// [[file:~/Workspace/Programming/gchemol/geometry.note::0e107ad6-6c63-45b9-9a8f-cb7cdd3e0777][0e107ad6-6c63-45b9-9a8f-cb7cdd3e0777]]
use super::*;

use nalgebra as na;

type Rotation3f = na::Rotation3<f64>;

mod qcprot;
mod quaternion;
// 0e107ad6-6c63-45b9-9a8f-cb7cdd3e0777 ends here

// [[file:~/Workspace/Programming/gchemol/geometry.note::ec9ee4ce-4967-4c41-bb0c-a823981b7631][ec9ee4ce-4967-4c41-bb0c-a823981b7631]]
/// The result of alignment defining how to superimpose.
#[derive(Clone, Debug)]
pub struct Superposition {
    /// superpostion rmsd
    pub rmsd: f64,
    /// translation vector
    pub translation: Vector3f,
    /// rotation matrix
    pub rotation_matrix: Matrix3f,
}

impl Superposition {
    /// Apply superposition to other structure
    pub fn apply(&self, conf: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut res = Vec::with_capacity(conf.len());
        for &v in conf {
            let v = Vector3f::from(v);
            let v = self.rotation_matrix * v + self.translation;
            res.push(v.into());
        }

        res
    }
}

/// Alignment of candidate structure onto the reference
#[derive(Clone, Debug)]
pub struct Alignment<'a> {
    /// The positions of the candidate structure
    positions: &'a [[f64; 3]],
}

impl<'a> Alignment<'a> {
    /// Construct from positions of the candidate to be aligned
    pub fn new(positions: &'a [[f64; 3]]) -> Self {
        Alignment {
            positions,
        }
    }

    /// Calculate Root-mean-square deviation of self with the reference coordinates
    ///
    /// Parameters
    /// ----------
    /// * reference: reference coordinates
    /// * weights  : weight of each point
    pub fn rmsd(&self, reference: &[[f64; 3]], weights: Option<&[f64]>) -> Result<f64> {
        // sanity check
        let npts = self.positions.len();
        if reference.len() != npts {
            bail!("points size mismatch!");
        }
        if weights.is_some() && weights.unwrap().len() != npts {
            bail!("weights size mismatch!");
        }

        // calculate rmsd
        // take the weight if any, or set it to 1.0
        let weights = weights.map_or_else(|| vec![1.0; npts].as_slice(), |ws| ws);
        let mut ws = 0.0f64;
        for i in 0..npts {
            let wi = weights[i];
            let dx = wi * (self.positions[i][0] - reference[i][0]);
            let dy = wi * (self.positions[i][1] - reference[i][1]);
            let dz = wi * (self.positions[i][2] - reference[i][2]);

            ws += dx.powi(2) + dy.powi(2) + dz.powi(2);
        }
        let ws = ws.sqrt();

        Ok(ws)
    }

    /// Superpose candidate structure onto reference structure which will be held fixed
    /// Return superposition struct
    ///
    /// Parameters
    /// ----------
    /// * reference: reference coordinates
    /// * weights  : weight of each point
    pub fn superpose(&mut self, reference: &[[f64; 3]], weights: Option<&[f64]>) -> Result<Superposition> {
        // calculate the RMSD & rotational matrix
        let (rmsd, trans, rot) = qcprot::calc_rmsd_rotational_matrix(&reference, &self.positions, weights);

        // return unit matrix if two structures are already close enough
        let rotation_matrix = if let Some(rot) = rot {
            Matrix3f::from_row_slice(&rot)
        } else {
            Matrix3f::identity()
        };

        // return superimposition result
        let sp = Superposition {
            rmsd,
            translation: trans.into(),
            rotation_matrix,
        };

        Ok(sp)
    }
}
// ec9ee4ce-4967-4c41-bb0c-a823981b7631 ends here

// [[file:~/Workspace/Programming/gchemol/geometry.note::128fc758-e836-41df-a94e-e90580bb73e3][128fc758-e836-41df-a94e-e90580bb73e3]]
#[test]
fn test_alignment() {
    // fragment a
    let (reference, candidate, weights) = qcprot::prepare_test_data();

    // construct alignment for superimposition
    let mut align = Alignment::new(&candidate);

    // alignment result
    let sp = align.superpose(&reference, Some(&weights)).unwrap();
    let rot = sp.rotation_matrix;

    // validation
    let rot_expected = Matrix3f::from_row_slice(&[
        0.77227551,    0.63510272,   -0.01533190,
        -0.44544846,    0.52413614,   -0.72584914,
        -0.45295276,    0.56738509,    0.68768304,
    ]);
    assert_relative_eq!(rot_expected, rot, epsilon=1e-4);
}

#[test]
fn test_alignment2() {
    use molecule::Molecule;

    // load test molecules
    let mol1 = Molecule::from_file("tests/files/alignment/reference.mol2").expect("alignment reference");
    let mol2 = Molecule::from_file("tests/files/alignment/candidate.mol2").expect("alignment candidate");

    // take the first 5 atoms for superposition
    let reference = &mol1.positions()[0..5];
    let candidate = &mol2.positions()[0..5];

    // align the candidate onto the reference
    let mut align = Alignment::new(&candidate);
    let sp = align.superpose(&reference, None).unwrap();

    // apply superposition to all atoms
    let new = sp.apply(&candidate);
    assert_relative_eq!(reference.to_dmatrix(), new.to_dmatrix(), epsilon=1e-3);
}
// 128fc758-e836-41df-a94e-e90580bb73e3 ends here
