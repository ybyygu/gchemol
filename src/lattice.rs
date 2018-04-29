// [[file:~/Workspace/Programming/gchemol/gchemol.note::891f59cf-3963-4dbe-a7d2-48279723b72e][891f59cf-3963-4dbe-a7d2-48279723b72e]]
//===============================================================================#
//   DESCRIPTION:  represents 3D periodic lattices
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-29 14:27>
//       UPDATED:  <2018-04-29 Sun 19:09>
//===============================================================================#

use nalgebra::{
    Vector3,               // A stack-allocated, 3-dimensional column vector.
    Matrix3,               // A stack-allocated, column-major, 3x3 square matrix
};

type Mat3D = Matrix3<f64>;
type Vec3D = Vector3<f64>;
// 891f59cf-3963-4dbe-a7d2-48279723b72e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b17e625d-352f-419e-9d10-a84fcdb9ff07][b17e625d-352f-419e-9d10-a84fcdb9ff07]]
#[derive(Debug, Clone)]
pub struct Lattice {
    /// internal translation matrix
    matrix: Mat3D,
    /// Lattice origin
    origin: Vec3D,
    /// Cached inverse of lattice matrix
    inv_matrix: Option<Mat3D>,
}

impl Default for Lattice {
    fn default() -> Self {
        Lattice {
            matrix     : Mat3D::identity(),
            origin     : Vec3D::zeros(),
            inv_matrix : None,
        }
    }
}

impl Lattice {
    /// using a cache to reduce the expensive matrix inversion calculations
    fn inv_matrix(&mut self) -> Mat3D {
        // make a readonly reference
        let matrix = self.matrix;
        let im = self.inv_matrix.get_or_insert_with(|| matrix.try_inverse().expect("bad matrix"));

        *im
    }

    pub fn new<T: Into<[[f64; 3]; 3]>>(tvs: T) -> Self {
        Lattice {
            matrix: Mat3D::from(tvs.into()),
            ..Default::default()
        }
    }

    /// Set cell origin in Cartesian coordinates
    pub fn set_origin(&mut self, loc: [f64; 3]) {
        self.origin = Vec3D::from(loc);
    }

    /// Lattice length parameters: a, b, c
    pub fn lengths(&self) -> (f64, f64, f64) {
        (
            self.matrix.column(0).norm(),
            self.matrix.column(1).norm(),
            self.matrix.column(2).norm()
        )
    }

    /// Lattice angle parameters in degrees
    pub fn angles(&self) -> (f64, f64, f64) {
        let va = self.matrix.column(0);
        let vb = self.matrix.column(1);
        let vc = self.matrix.column(2);
        (
            vb.angle(&vc).to_degrees(),
            va.angle(&vc).to_degrees(),
            va.angle(&vb).to_degrees(),
        )
    }

    /// Get cell origin in Cartesian coordinates
    pub fn origin(&self) -> [f64; 3] {
        self.origin.into()
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    pub fn to_frac(&mut self, p: [f64; 3]) -> [f64; 3] {
        let im = self.inv_matrix();
        let v = Vec3D::from(p);
        let fs = im*v;
        fs.into()
    }

    /// Returns the cartesian coordinates given fractional coordinates.
    pub fn to_cart(&self, p: [f64; 3]) -> [f64; 3] {
        let v = Vec3D::from(p);
        let fs = self.matrix*v + self.origin;

        fs.into()
    }
}
// b17e625d-352f-419e-9d10-a84fcdb9ff07 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::4bc21235-f285-4976-a32a-b33506381b58][4bc21235-f285-4976-a32a-b33506381b58]]
#[test]
fn test_lattice_init() {
    let mut lat = Lattice::default();
    let loc = [1.0, 2.0, 3.0];
    lat.set_origin(loc);
    assert_eq!(loc, lat.origin());

    let lat = Lattice::new([[ 15.3643,   0.    ,   0.    ],
                            [  4.5807,  15.5026,   0.    ],
                            [  0.    ,   0.    ,  17.4858]]);

    let (a, b, c) = lat.lengths();

    assert_relative_eq!(a, 15.3643, epsilon=1e-4);
    assert_relative_eq!(b, 16.1652, epsilon=1e-4);
    assert_relative_eq!(c, 17.4858, epsilon=1e-4);

    let (alpha, beta, gamma) = lat.angles();
    assert_relative_eq!(alpha, 90.0, epsilon=1e-4);
    assert_relative_eq!(beta, 90.0, epsilon=1e-4);
    assert_relative_eq!(gamma, 73.5386, epsilon=1e-4);
}

#[test]
fn test_lattice() {
    // ovito/tests/files/LAMMPS/multi_sequence_1.dump
    let mut lat = Lattice::new([[5.09, 0.00, 0.00],
                                [0.00, 6.74, 0.00],
                                [0.00, 0.00, 4.53]]);

    let fs = lat.to_frac([2.1832, 1.6850, 3.8505]);
    assert_relative_eq!(fs[0], 0.4289, epsilon=1e-3);
    assert_relative_eq!(fs[1], 0.2500, epsilon=1e-3);
    assert_relative_eq!(fs[2], 0.8500, epsilon=1e-3);
    let fs = lat.to_frac([6.9068, 5.0550, 0.6795]);
    assert_relative_eq!(fs[0], 1.3569, epsilon=1e-3);
    let fs = lat.to_frac([4.3618, 5.0550, 1.5855]);
    assert_relative_eq!(fs[2], 0.3500, epsilon=1e-3);

    let coords = lat.to_cart([0.4289, 0.2500, 0.8500]);
    assert_relative_eq!(coords[0], 2.1832, epsilon=1e-3);
    assert_relative_eq!(coords[1], 1.6850, epsilon=1e-3);
    assert_relative_eq!(coords[2], 3.8505, epsilon=1e-3);
}
// 4bc21235-f285-4976-a32a-b33506381b58 ends here
