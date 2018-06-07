// [[file:~/Workspace/Programming/gchemol/gchemol.note::891f59cf-3963-4dbe-a7d2-48279723b72e][891f59cf-3963-4dbe-a7d2-48279723b72e]]
//===============================================================================#
//   DESCRIPTION:  Represents 3D periodic lattices
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  ---
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-29 14:27>
//       UPDATED:  <2018-06-07 Thu 17:20>
//===============================================================================#

use nalgebra::{
    Vector3,               // A stack-allocated, 3-dimensional column vector.
    Matrix3,               // A stack-allocated, column-major, 3x3 square matrix
};

type Mat3D = Matrix3<f64>;
type Vec3D = Vector3<f64>;
// 891f59cf-3963-4dbe-a7d2-48279723b72e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b17e625d-352f-419e-9d10-a84fcdb9ff07][b17e625d-352f-419e-9d10-a84fcdb9ff07]]
/// Periodic 3D lattice
#[derive(Debug, Clone, Copy)]
pub struct Lattice {
    /// internal translation matrix
    matrix: Mat3D,
    /// Lattice origin
    origin: Vec3D,

    /// Cached inverse of lattice matrix
    inv_matrix: Option<Mat3D>,

    /// Cached volume of the unit cell.
    volume: Option<f64>,

    /// The perpendicular widths of the unit cell on each direction,
    /// i.e. the distance between opposite faces of the unit cell
    widths: Option<[f64; 3]>,

    /// Cached cell lengths parameters
    lengths: Option<[f64; 3]>,

    /// Cached cell angles parameters
    angles: Option<[f64; 3]>,
}

// matrix inversion
fn get_inv_matrix(matrix: Mat3D) -> Mat3D {
    matrix.try_inverse().expect("bad matrix")
}

// cell volume
fn get_cell_volume(mat: Mat3D) -> f64 {
    let va = mat.column(0);
    let vb = mat.column(1);
    let vc = mat.column(2);
    va.dot(&vb.cross(&vc))
}

// return cell length parameters
fn get_cell_lengths(mat: Mat3D) -> [f64; 3] {
    [
        mat.column(0).norm(),
        mat.column(1).norm(),
        mat.column(2).norm()
    ]
}

// return cell angle parameters in degrees
fn get_cell_angles(mat: Mat3D) -> [f64; 3] {
    let va = mat.column(0);
    let vb = mat.column(1);
    let vc = mat.column(2);
    [
        vb.angle(&vc).to_degrees(),
        va.angle(&vc).to_degrees(),
        va.angle(&vb).to_degrees(),
    ]
}

impl Default for Lattice {
    fn default() -> Self {
        Lattice {
            matrix: Mat3D::identity(),
            origin: Vec3D::zeros(),

            inv_matrix: None,
            volume: None,
            widths: None,
            lengths: None,
            angles: None,
        }
    }
}

impl Lattice {
    pub fn new<T: Into<[[f64; 3]; 3]>>(tvs: T) -> Self {
        Lattice {
            matrix: Mat3D::from(tvs.into()),
            ..Default::default()
        }
    }

    /// using a cache to reduce the expensive matrix inversion calculations
    fn inv_matrix(&mut self) -> Mat3D {
        // make a readonly reference
        let matrix = self.matrix;
        let im = self.inv_matrix.get_or_insert_with(|| matrix.try_inverse().expect("bad matrix"));

        *im
    }

    fn get_cell_widths(&mut self)  -> [f64; 3] {
        let volume = self.volume();
        let [van, vbn, vcn] = self.lengths();

        let wa = volume / (vbn*vcn);
        let wb = volume / (vcn*van);
        let wc = volume / (van*vbn);

        [wa, wb, wc]
    }

    pub fn widths(&mut self) -> [f64; 3] {
        if let Some(ws) = self.widths {
            return ws;
        } else {
            let ws = self.get_cell_widths();
            self.widths = Some(ws);

            ws
        }
    }

    /// Return the volume of the unit cell
    /// the cache will be updated if necessary
    pub fn volume(&mut self) -> f64 {
        // make a read-only reference
        let mat = self.matrix;
        let volume = self.volume.get_or_insert_with(|| get_cell_volume(mat));

        *volume
    }

    /// Construct lattice from lattice parameters
    /// Unit cell angles in degrees, lengths in Angstrom
    pub fn from_params(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let alpha = alpha.to_radians();
        let beta  = beta.to_radians();
        let gamma = gamma.to_radians();

        let acos = alpha.cos();
        let bcos = beta.cos();
        let gcos = gamma.cos();
        let gsin = gamma.sin();
        let v = (1.
                 - acos.powi(2)
                 - bcos.powi(2)
                 - gcos.powi(2)
                 + 2.0 * acos * bcos * gcos).sqrt();

        let va = [a,
                  0.0,
                  0.0];

        let vb = [b*gcos,
                  b*gsin,
                  0.0
                  ];

        let vc = [c*bcos,
                  c*(acos - bcos*gcos)/gsin,
                  c*v/gsin];

        Lattice::new([va, vb, vc])
    }

    /// Set cell origin in Cartesian coordinates
    pub fn set_origin(&mut self, loc: [f64; 3]) {
        self.origin = Vec3D::from(loc);
    }

    /// Lattice length parameters: a, b, c
    pub fn lengths(&mut self) -> [f64; 3] {
        let mat = self.matrix;
        let lengths = self.lengths.get_or_insert_with(|| get_cell_lengths(mat));

        [
            lengths[0],
            lengths[1],
            lengths[2],
        ]
    }

    /// Lattice angle parameters in degrees
    pub fn angles(&mut self) -> [f64; 3] {
        let mat = self.matrix;
        let angles = self.angles.get_or_insert_with(|| get_cell_angles(mat));

        [
            angles[0],
            angles[1],
            angles[2],
        ]
    }

    // FIXME: cell widths
    /// Scale Lattice by a positive constant
    pub fn scale_by(&mut self, v: f64) {
        debug_assert!(v > 0.);
        self.matrix *= v;

        // reset caches
        self.inv_matrix = None;
        self.volume = None;
        self.widths = None;
        self.lengths = None;
        self.angles = None;
    }

    /// Get cell origin in Cartesian coordinates
    pub fn origin(&self) -> [f64; 3] {
        self.origin.into()
    }

    /// Returns the fractional coordinates given cartesian coordinates.
    pub fn to_frac(&mut self, p: [f64; 3]) -> [f64; 3] {
        let im = self.inv_matrix();
        let v = Vec3D::from(p);
        let fs = im * (v - self.origin);
        fs.into()
    }

    /// Returns the cartesian coordinates given fractional coordinates.
    pub fn to_cart(&self, p: [f64; 3]) -> [f64; 3] {
        let v = Vec3D::from(p);
        let fs = self.matrix*v + self.origin;

        fs.into()
    }

    /// Lattice vector a
    pub fn vector_a(&self) -> [f64; 3] {
        self.matrix.column(0).transpose().into()
    }

    /// Lattice vector b
    pub fn vector_b(&self) -> [f64; 3] {
        self.matrix.column(1).transpose().into()
    }

    /// Lattice vector c
    pub fn vector_c(&self) -> [f64; 3] {
        self.matrix.column(2).transpose().into()
    }

    /// Lattice vectors
    pub fn vectors(&self) -> [[f64; 3]; 3] {
        self.matrix.into()
    }

    /// Check if lattice is orthorhombic
    pub fn is_orthorhombic(&self) -> bool {
        let diag = self.matrix.diagonal();
        let m = Matrix3::from_diagonal(&diag);
        m == self.matrix
    }
}
// b17e625d-352f-419e-9d10-a84fcdb9ff07 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::83cff231-cc63-4077-b07e-a26a2c2b906d][83cff231-cc63-4077-b07e-a26a2c2b906d]]
use std::f64;

impl Lattice {
    /// Return the relevant periodic images required for neighborhood search
    /// within cutoff radius
    pub fn relevant_images(&mut self, radius: f64) -> Vec<Vector3<f64>> {
        let ns = self.n_min_images(radius);
        let na = ns[0] as isize;
        let nb = ns[1] as isize;
        let nc = ns[2] as isize;

        let mut images = vec![];
        for i in (-na)..(na+1) {
            for j in (-nb)..(nb+1) {
                for k in (-nc)..(nc+1) {
                    let v = Vec3D::from([i as f64, j as f64, k as f64]);
                    images.push(v);
                }
            }
        }

        images
    }

    /// Return the minimal number of images for neighborhood search on each cell
    /// direction within cutoff radius
    fn n_min_images(&mut self, radius: f64) -> [usize; 3]{
        let mut ns = [0; 3];

        for (i, &w) in self.widths().iter().enumerate() {
            let n = (radius / w).ceil();
            ns[i] = n as usize;
        }

        ns
    }

    /// Return the distance between two points computed using the minimum image
    /// convention.
    ///
    /// Reference
    /// ---------
    /// (1) Tuckerman, M. E. Statistical Mechanics: Theory and Molecular
    /// Simulation, 1 edition.; Oxford University Press: Oxford ; New York,
    /// 2010.
    fn distance_tuckerman(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        let pij = [pj[0] - pi[0],
                   pj[1] - pi[1],
                   pj[2] - pi[2]];

        // apply minimum image convention on the scaled coordinates
        let mut fij = self.to_frac(pij);
        for x in 0..3 {
            fij[x] -= fij[x].round();
        }

        // transform back to cartesian coordinates
        let pij = self.to_cart(fij);
        Vec3D::from(pij).norm()
    }


    /// Return the shortest distance between `pi` and the periodic images of `pj`
    /// This algorithm will loop over all relevant images
    fn distance_brute_force(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        // the cutoff radius for finding relevant images
        let cutoff = self.distance_tuckerman(pi, pj);
        let relevant_images = self.relevant_images(cutoff);

        let mut distance = f64::MAX;
        for v in relevant_images {
            let dd = self.to_cart(v.into());
            let ipij = [
                pj[0] + dd[0] - pi[0],
                pj[1] + dd[1] - pi[1],
                pj[2] + dd[2] - pi[2],
            ];

            let d = Vec3D::from(ipij).norm();
            if d < distance {
                distance = d;
            }
        }

        distance
    }

    /// Return the shortest distance between `pi` and the periodic images of `pj`
    pub fn distance(&mut self, pi: [f64; 3], pj: [f64; 3]) -> f64 {
        if self.is_orthorhombic() {
            self.distance_tuckerman(pi, pj)
        } else {
            self.distance_brute_force(pi, pj)
        }
    }
}
// 83cff231-cc63-4077-b07e-a26a2c2b906d ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::4bc21235-f285-4976-a32a-b33506381b58][4bc21235-f285-4976-a32a-b33506381b58]]
#[test]
fn test_lattice_construct() {
    let mut lat = Lattice::default();
    let loc = [1.0, 2.0, 3.0];
    lat.set_origin(loc);
    assert_eq!(loc, lat.origin());

    let mut lat = Lattice::new([[ 15.3643,   0.    ,   0.    ],
                                [  4.5807,  15.5026,   0.    ],
                                [  0.    ,   0.    ,  17.4858]]);

    let [a, b, c] = lat.lengths();
    assert_eq!(false, lat.is_orthorhombic());

    assert_relative_eq!(a, 15.3643, epsilon=1e-4);
    assert_relative_eq!(b, 16.1652, epsilon=1e-4);
    assert_relative_eq!(c, 17.4858, epsilon=1e-4);

    let [alpha, beta, gamma] = lat.angles();
    assert_relative_eq!(alpha, 90.0, epsilon=1e-4);
    assert_relative_eq!(beta, 90.0, epsilon=1e-4);
    assert_relative_eq!(gamma, 73.5386, epsilon=1e-4);

    let mut lat = Lattice::from_params(a, b, c, alpha, beta, gamma);
    assert_eq!([a, b, c], lat.lengths());
    assert_eq!([alpha, beta, gamma], lat.angles());
}

#[test]
fn test_lattice_neighborhood() {
    let mut lat = Lattice::new([
        [ 18.256,   0.   ,   0.   ],
        [  0.   ,  20.534,   0.   ],
        [  0.   ,   0.   ,  15.084],
    ]);
    assert_eq!(true, lat.is_orthorhombic());

    assert_eq!([1, 1, 1], lat.n_min_images(9.));
    assert_eq!([2, 1, 2], lat.n_min_images(19.));
    assert_eq!([2, 1, 2], lat.n_min_images(20.));
    assert_eq!([2, 2, 2], lat.n_min_images(20.6));

    let expected = [
        Vec3D::new(-1.0, -1.0, -1.0),
        Vec3D::new(-1.0, -1.0,  0.0),
        Vec3D::new(-1.0, -1.0,  1.0),
        Vec3D::new(-1.0,  0.0, -1.0),
        Vec3D::new(-1.0,  0.0,  0.0),
        Vec3D::new(-1.0,  0.0,  1.0),
        Vec3D::new(-1.0,  1.0, -1.0),
        Vec3D::new(-1.0,  1.0,  0.0),
        Vec3D::new(-1.0,  1.0,  1.0),
        Vec3D::new( 0.0, -1.0, -1.0),
        Vec3D::new( 0.0, -1.0,  0.0),
        Vec3D::new( 0.0, -1.0,  1.0),
        Vec3D::new( 0.0,  0.0, -1.0),
        Vec3D::new( 0.0,  0.0,  0.0),
        Vec3D::new( 0.0,  0.0,  1.0),
        Vec3D::new( 0.0,  1.0, -1.0),
        Vec3D::new( 0.0,  1.0,  0.0),
        Vec3D::new( 0.0,  1.0,  1.0),
        Vec3D::new( 1.0, -1.0, -1.0),
        Vec3D::new( 1.0, -1.0,  0.0),
        Vec3D::new( 1.0, -1.0,  1.0),
        Vec3D::new( 1.0,  0.0, -1.0),
        Vec3D::new( 1.0,  0.0,  0.0),
        Vec3D::new( 1.0,  0.0,  1.0),
        Vec3D::new( 1.0,  1.0, -1.0),
        Vec3D::new( 1.0,  1.0,  0.0),
        Vec3D::new( 1.0,  1.0,  1.0)];

    let images = lat.relevant_images(3.0);
    assert_eq!(expected.len(), images.len());
    assert_eq!(expected[1][2], images[1][2]);
}

#[test]
fn test_lattice_volume() {
    let vts = [
        [ 5.,  0.,  0.],
        [ 5.,  5.,  0.],
        [ 1.,  0.,  5.]
    ];

    let mut lat = Lattice::new(vts);
    assert_eq!(vts, lat.vectors());
    assert_eq!(vts[0], lat.vector_a());
    assert_eq!(vts[1], lat.vector_b());
    assert_eq!(vts[2], lat.vector_c());

    assert_relative_eq!(125.0, lat.volume(), epsilon=1e-4);
    lat.scale_by(4.);
    assert_relative_eq!(8000.0, lat.volume(), epsilon=1e-4);
}

#[test]
fn test_lattice_frac_cart() {
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

#[test]
fn test_lattice_mic_distance() {
    let mut lat = Lattice::new([
        [5.0, 0.0, 0.0],
        [1.0, 5.0, 0.0],
        [1.0, 1.0, 5.0],
    ]);

    // the shortest distance: 2.61383
    let d = lat.distance_tuckerman([0.; 3], [-0.94112, -4.34823, 2.53058]);
    assert_relative_eq!(2.66552, d, epsilon=1e-4);
    let d = lat.distance_brute_force([0.; 3], [-0.94112, -4.34823, 2.53058]);
    assert_relative_eq!(2.61383, d, epsilon=1e-4);

    // the shortest distance: 2.53575
    let d = lat.distance_tuckerman([0.; 3], [-2.46763, 0.57717, 0.08775]);
    assert_relative_eq!(2.59879, d, epsilon=1e-4);
    let d = lat.distance_brute_force([0.; 3], [-2.46763, 0.57717, 0.08775]);
    assert_relative_eq!(2.53575, d, epsilon=1e-4);
}
// 4bc21235-f285-4976-a32a-b33506381b58 ends here
