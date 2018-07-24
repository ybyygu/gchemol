// [[file:~/Workspace/Programming/gchemol/geometry.note::f2ed012d-0e00-4288-b59e-0cb61f7921c2][f2ed012d-0e00-4288-b59e-0cb61f7921c2]]
/// Providing simple statistics methods (min, max, mean, var, ...) for [f64]
pub use test::stats::Stats;
extern crate nalgebra;
use nalgebra as na;

/// Vector in 3D space
pub type Vector3f = na::Vector3<f64>;

/// 3xN matrix storing a list of 3D vectors
pub type Vector3fVec = na::Matrix<f64, na::U3, na::Dynamic, na::MatrixVec<f64, na::U3, na::Dynamic>>;

/// Cartesian 3D coordinates
pub type Position = Vector3f;
/// A list of 3D position vectors
pub type Positions = Vector3fVec;

/// Dynamic array containing float numbers
pub type FloatVec = na::DVector<f64>;

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

    fn to_dmatrix(&self) -> Vector3fVec {
        // FIXME: performance
        let r = self.ravel();
        Vector3fVec::from_column_slice(self.len(), &r)
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

mod transform;
mod random;

pub use self::transform::*;
pub use self::random::*;
// f2ed012d-0e00-4288-b59e-0cb61f7921c2 ends here
