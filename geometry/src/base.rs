// types

// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::*types][types:1]]
use quicli::prelude::*;

/// Providing simple statistics methods (min, max, mean, var, ...) for [f64]
pub use crate::test::stats::Stats;
use nalgebra as na;

/// Vector in 3D space
pub type Vector3f = na::Vector3<f64>;
pub type Matrix3f = na::Matrix3<f64>;
pub type Matrix4f = na::Matrix4<f64>;
pub type DMatrixf = na::DMatrix<f64>;

/// 3xN matrix storing a list of 3D vectors
pub type Vector3fVec = na::Matrix<f64, na::U3, na::Dynamic, na::MatrixVec<f64, na::U3, na::Dynamic>>;

/// Cartesian 3D coordinates
pub type Position = Vector3f;
/// A list of 3D position vectors
pub type Positions = Vector3fVec;

/// Dynamic array containing float numbers
pub type FloatVec = na::DVector<f64>;
// types:1 ends here

// traits

// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::*traits][traits:1]]
/// A trait provides useful tools for Vec<f64> type.
pub trait VecFloatMath {
    /// Convert to nalgebra dynamic 1xN vector
    fn to_dvector(&self) -> FloatVec;

    /// 1xN Vector norm
    fn norm(&self) -> f64;
}

impl VecFloatMath for [f64] {
    fn to_dvector(&self) -> FloatVec {
        FloatVec::from_column_slice(self.len(), &self)
    }

    fn norm(&self) -> f64 {
        let mut s = 0.0;
        for &x in self {
            s += x.powi(2);
        }

        s.sqrt()
    }
}

/// A trait provides useful tools for Vec<[f64; 3]> type.
pub trait VecFloat3Math {
    /// Return the Frobenius norm (matrix norm) defined as the square root of
    /// the sum of the absolute squares of its elements.
    fn norm(&self) -> f64;

    /// return the norms of a list of 3D vectors
    fn norms(&self) -> Vec<f64>;

    /// Return a 1-D array, containing the elements of 3xN array
    fn ravel(&self) -> Vec<f64>;

    /// Convert to 3xN dynamic matrix
    fn to_dmatrix(&self) -> Vector3fVec;

    /// Convert to Vec of Vector3f
    fn to_vectors(&self) -> Vec<Vector3f>;

    /// Return the center of geometry (COG)
    fn center_of_geometry(&self) -> Position;

    /// Return the center of geometry
    ///
    /// https://en.wikipedia.org/wiki/Centroid
    fn centroid(&self) -> Position {
        self.center_of_geometry()
    }

    /// Return weighted center of geometry (COM)
    fn center_of_mass(&self, masses: &[f64]) -> Result<Position>;

    /// Return mirror inverted structure
    fn mirror_inverted(&self) -> Positions;

    /// Return point inverted structure
    fn point_inverted(&self) -> Positions;

    /// Return distance matrix
    fn distance_matrix(&self) -> DMatrixf;
}

// impl VecFloat3Math for Vec<[f64; 3]> {
impl VecFloat3Math for [[f64; 3]] {
    fn ravel(&self) -> Vec<f64> {
        let n = self.len();
        let mut r = Vec::with_capacity(3 * n);

        for i in 0..n {
            for j in 0..3 {
                r.push(self[i][j]);
            }
        }

        r
    }

    fn norm(&self) -> f64 {
        let mut d2: f64 = 0.0;
        for i in 0..self.len() {
            for j in 0..3 {
                d2 += self[i][j].powi(2);
            }
        }
        d2.sqrt()
    }

    fn norms(&self) -> Vec<f64> {
        let n = self.len();
        let mut norms = Vec::with_capacity(n);

        for i in 0..n {
            let mut l = 0.0;
            for j in 0..3 {
                let vij = self[i][j];
                l += vij.powi(2);
            }

            norms.push(l.sqrt());
        }

        norms
    }

    fn to_vectors(&self) -> Vec<Vector3f> {
        self.iter().map(|&a| Vector3f::from(a)).collect()
    }

    fn to_dmatrix(&self) -> Vector3fVec {
        // FIXME: performance
        let r = self.ravel();
        Vector3fVec::from_column_slice(self.len(), &r)
    }

    fn center_of_mass(&self, masses: &[f64]) -> Result<Position> {
        weighted_center_of_geometry(&self, &masses)
    }

    fn center_of_geometry(&self) -> Position {
        let n = self.len();
        let weights: Vec<_> = (0..n).map(|_ |1.0).collect();
        weighted_center_of_geometry(&self, &weights).expect("center of geometry")
    }

    fn mirror_inverted(&self) -> Positions {
        let m = self.to_dmatrix();
        let r = na::Matrix3::from_diagonal(&[1.0, 1.0, -1.0].into());
        r * m
    }

    fn point_inverted(&self) -> Positions {
        let m = self.to_dmatrix();
        let r = na::Matrix3::from_diagonal(&[-1.0, -1.0, -1.0].into());
        r * m
    }

    fn distance_matrix(&self) -> DMatrixf {
        let n = self.len();

        let mut distances = DMatrixf::zeros(n, n);
        for i in 0..n {
            for j in 0..i{
                let d = euclidean_distance(self[i], self[j]);
                distances[(i, j)] = d;
                distances[(j, i)] = d;
            }
        }

        distances
    }
}

#[test]
fn test_vec_math() {
    let a = vec![1.0, 2.0, 3.0];
    let x = a.to_dvector();
    assert_relative_eq!(x.norm(), a.norm(), epsilon=1e-3);

    let positions = [[-0.131944, -0.282942,  0.315957],
                     [ 0.40122 , -1.210646,  0.315957],
                     [-1.201944, -0.282942,  0.315957],
                     [ 0.543331,  0.892036,  0.315957],
                     [ 0.010167,  1.819741,  0.315957],
                     [ 1.613331,  0.892036,  0.315957]];

    let n = positions.norm();
    let m = positions.to_dmatrix();
    assert_relative_eq!(n, m.norm(), epsilon=1e-4);
    let n = positions.to_vec().norm();

    let x = positions.ravel();
    assert_eq!(positions.len() * 3, x.len());

    let x = positions.norms().max();
    assert_relative_eq!(1.8704, x, epsilon=1e-4);
}
// traits:1 ends here

// functions

// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::*functions][functions:1]]
#[inline]
pub fn euclidean_distance(p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let mut d2 = 0.0;
    for v in 0..3 {
        let dv = p2[v] - p1[v];
        d2 += dv*dv;
    }

    d2.sqrt()
}

/// Return the geometric center
pub fn weighted_center_of_geometry(positions: &[[f64; 3]], weights: &[f64]) -> Result<Vector3f> {
    let npts = positions.len();
    let mut pc = [0.0; 3];

    // sanity check
    if npts != weights.len() {
        bail!("size inconsistent!");
    }

    // deviding by zero?
    let mut wsum = weights.sum();
    if wsum < 1e-6 {
        error!("weird weight sum: {:?}", wsum);
    }

    for i in 0..npts {
        for j in 0..3 {
            pc[j] += weights[i] * positions[i][j];
        }
    }

    for i in 0..3 {
        pc[i] /= wsum;
    }

    Ok(Vector3f::from(pc))
}

#[test]
fn test_weighted_center_of_geometry() {
    // points
    let frag = vec![
        [ -2.803, -15.373, 24.556],
        [  0.893, -16.062, 25.147],
        [  1.368, -12.371, 25.885],
        [ -1.651, -12.153, 28.177],
        [ -0.440, -15.218, 30.068],
        [  2.551, -13.273, 31.372],
        [  0.105, -11.330, 33.567],
    ];

    // weights
    let natoms = frag.len();
    let masses: Vec<_> = (0..natoms).map(|v| v as f64 + 1.0).collect();

    // expected results
    let expected = Vector3f::new(0.3687142857142857, -13.15214285714286, 29.955499999999997);
    let pc = frag.center_of_mass(&masses).expect("geometry: com");
    relative_eq!(pc, expected, epsilon=1e-6);
}
// functions:1 ends here
