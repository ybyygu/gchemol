// [[file:~/Workspace/Programming/gchemol/geometry.note::0e107ad6-6c63-45b9-9a8f-cb7cdd3e0777][0e107ad6-6c63-45b9-9a8f-cb7cdd3e0777]]
use quicli::prelude::*;

use nalgebra as na;
use super::*;

mod qcprot;

/// The result of alignment defining how to superimpose.
#[derive(Clone, Debug)]
pub struct Superposition {
    // superpostion rmsd
    rmsd: f64,
    // translation vector
    translation: Vector3f,
    // rotation matrix
    rotation_matrix: Option<Matrix3f>,
}

/// Alignment of a set of points onto the reference points
#[derive(Clone, Debug)]
pub struct Alignment {
    /// The position vector of each point
    positions: Vec<[f64; 3]>,

    /// The weight of each point
    weights: Option<FloatVec>,
}

impl Alignment {
    /// Construct from positions of the points to be aligned
    pub fn new(positions: &[[f64; 3]]) -> Self {
        let positions = positions.into();

        Alignment {
            positions,
            weights: None,
        }
    }

    /// Superpose current points onto reference points which will be held fixed
    /// Return superposition RMSD
    /// Parameters
    /// ----------
    /// * reference: reference coordinates
    /// * weights  : weight of each point
    pub fn superpose(&self, reference: &[[f64; 3]], weights: Option<&[f64]>) -> Result<Superposition> {
        // calculate the RMSD & rotational matrix
        let (rmsd, trans, rot9) = qcprot::calc_rmsd_rotational_matrix(&reference, &self.positions, weights);

        // return superimposition result
        let sp = Superposition {
            rmsd,
            translation: trans.into(),
            rotation_matrix: rot9.map(|a| Matrix3f::from_column_slice(&a)),
        };

        Ok(sp)
    }
}

#[test]
fn test_alignment() {
    // fragment a
    let (reference, candidate, weights) = qcprot::prepare_test_data();

    let rot_expected = Matrix3f::from_column_slice(&[
        0.77227551,    0.63510272,   -0.01533190,
        -0.44544846,    0.52413614,   -0.72584914,
        -0.45295276,    0.56738509,    0.68768304,
    ]);

    let mut align = Alignment::new(&candidate);

    // alignment result
    let sp = align.superpose(&reference, Some(&weights)).unwrap();
    let rot = sp.rotation_matrix.expect("superposition: rotation matrix");
    assert_relative_eq!(rot_expected, rot, epsilon=1e-4);
}
// 0e107ad6-6c63-45b9-9a8f-cb7cdd3e0777 ends here
