// [[file:~/Workspace/Programming/gchemol/gchemol.note::f2ed012d-0e00-4288-b59e-0cb61f7921c2][f2ed012d-0e00-4288-b59e-0cb61f7921c2]]
use Point3D;
use Points;

use nalgebra::{
    Vector3,
    Rotation3,
};

pub trait VectorMath {
    /// Return the norm of 1D vector
    fn norm(&self) -> f64;
}

impl VectorMath for [f64] {
    fn norm(&self) -> f64 {
        let mut s = 0.0;
        for &x in self {
            s += x.powi(2);
        }

        s.sqrt()
    }
}

#[test]
fn test_vector_math() {
    let x = vec![0.1, 0.1, 0.1];
    let x = x.norm();
    assert_relative_eq!(0.1732, x, epsilon=1e-4);
}

pub trait PointMath {
    /// Return 3D vector norm.
    fn norm(&self) -> f64;

    /// Return the distance between two 3D points
    fn distance(&self, other: Point3D) -> f64;

    /// Return the dot product of two 3D vectors.
    fn vdot(&self, other: Point3D) -> f64;

    /// Return the length of vector.
    fn length(&self) -> f64 {
        self.norm()
    }

    fn x(&self) -> f64;
    fn y(&self) -> f64;
    fn z(&self) -> f64;
}

impl PointMath for Point3D {
    fn x(&self) -> f64 {
        self[0]
    }

    fn y(&self) -> f64 {
        self[1]
    }

    fn z(&self) -> f64 {
        self[2]
    }

    fn norm(&self) -> f64 {
        let mut d2: f64 = 0.0;
        for v in 0..3 {
            d2 += self[v].powi(2);
        }
        d2.sqrt()
    }

    fn distance(&self, other: Point3D) -> f64 {
        euclidean_distance(*self, other)
    }

    /// Return the dot product of two 3D vectors.
    fn vdot(&self, other: Point3D) -> f64 {
        let mut s = 0.0;
        for i in 0..self.len() {
            s += self[i] * other[i];
        }

        s
    }
}

pub trait PointsMath {
    fn vdot(&self, other: &Points) -> f64;

    /// Return the Frobenius norm (matrix norm) defined as the square root of
    /// the sum of the absolute squares of its elements.
    fn norm(&self) -> f64;

    /// return the norms of a list of 3D vectors
    fn norms(&self) -> Vec<f64>;

    /// Return the scaled positions by a constant
    fn scaled(&self, f: f64) -> Vec<Point3D>;
}

impl PointsMath for Points {
    fn vdot(&self, other: &Points) -> f64 {
        let n = self.len();
        debug_assert!(n == other.len());

        let mut vret = 0.0;
        for i in 0..n {
            for j in 0..3 {
                let vij1 = self[i][j];
                let vij2 = other[i][j];
                vret += vij1 * vij2;
            }
        }

        vret
    }

    fn scaled(&self, f: f64) -> Vec<Point3D> {
        let n = self.len();
        let mut pts = Vec::with_capacity(n);
        for i in 0..n {
            let mut p = [0.0; 3];
            for j in 0..3 {
                p[j] = self[i][j] * f;
            }
            pts.push(p);
        }
        pts
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
}

#[test]
fn test_geometry_trait() {
    let p1 = [0.0, 0.1, 0.0];
    let p2 = [0.0, 0.1, 0.0];
    let x = p1.norm();
    assert_eq!(x, p1.length());

    assert_eq!(0.1, p1.y());

    let d = p1.distance(p2);
    assert_relative_eq!(0.0, d, epsilon=1e-4);

    let v1 = [1.0, 2.0, 3.0];
    let v2 = [3.0, 4.0, 5.0];
    let d = v1.vdot(v2);
    assert_relative_eq!(26.0, d, epsilon=1e-4);

    let pts1 = vec![[1.0, 2.0, 3.0], [1.0, 1.0, 1.0]];
    let pts2 = vec![[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]];
    let d = pts1.vdot(&pts2);
    assert_relative_eq!(29.0, d, epsilon=1e-4);

    let n = pts1.norm();
    assert_relative_eq!(4.123126, n, epsilon=1e-4);
}

#[inline]
pub fn euclidean_distance(p1: Point3D, p2: Point3D) -> f64 {
    let mut d2 = 0.0;
    for v in 0..3 {
        let dv = p2[v] - p1[v];
        d2 += dv*dv;
    }

    d2.sqrt()
}

/// check if any pair of points come too close
pub fn close_contact(points: &Points) -> bool {
    let cutoff = 0.4;

    let npts = points.len();
    for i in 0..npts {
        for j in (i+1)..npts {
            let p1 = points[i];
            let p2 = points[j];
            let dx = p2[0] - p1[0];
            let dy = p2[1] - p1[1];
            let dz = p2[2] - p1[2];
            let d2 = dx*dx + dy*dy + dz*dz;
            if d2 <= cutoff {
                return true
            }
        }
    }

    false
}

/// Return all distances between any pair of points
pub fn get_distance_matrix(points: Points) -> Vec<Vec<f64>>{
    let npts = points.len();

    // fill distance matrix
    let mut distmat = vec![];
    for i in 0..npts {
        let mut dijs = vec![];
        for j in 0..npts {
            let dij = euclidean_distance(points[i], points[j]);
            dijs.push(dij);
        }
        distmat.push(dijs);
    }

    distmat
}
// f2ed012d-0e00-4288-b59e-0cb61f7921c2 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::26f9be7f-1dbd-4ac0-9cdf-4759ede5d338][26f9be7f-1dbd-4ac0-9cdf-4759ede5d338]]
/// Translate all points to a new location
pub fn translate(points: &mut Points, loc: Point3D) {
    for i in 0..points.len() {
        for v in 0..3 {
            points[i][v] += loc[v];
        }
    }
}
// 26f9be7f-1dbd-4ac0-9cdf-4759ede5d338 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c7248d32-74e1-47a8-8a78-56fef9ca529c][c7248d32-74e1-47a8-8a78-56fef9ca529c]]
/// rotate coordinates about x axis in radian
pub fn rotate_about_x_axis(points: &Points, angle: f64, center: Point3D) -> Points {
    let axis = Vector3::x_axis();
    let r = Rotation3::from_axis_angle(&axis, angle);

    let mut rpoints = vec![];
    let center = Vector3::from(center);
    for &p in points.iter() {
        let v = Vector3::from(p) - center;
        let t: Point3D = (r*v + center).into();
        rpoints.push(t);
    }

    rpoints
}
// c7248d32-74e1-47a8-8a78-56fef9ca529c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::9413e1bc-8f8f-4b07-b305-9d2911afabc6][9413e1bc-8f8f-4b07-b305-9d2911afabc6]]
use rand::{self, Rng};
use rand::distributions::{Range, Normal};

/// create a random point within a sphere
/// References
/// ----------
/// https://stackoverflow.com/a/5408344
pub fn rand_point_within_sphere(radius: f64) -> Point3D {
    debug_assert!(radius.is_sign_positive(), "sphere radius cannot be negative: {:?}", radius);
    let mut rng = rand::thread_rng();
    let range = Range::new(-radius, radius);

    // using the discarding method, which is simple and also fast
    let radius2 = radius*radius;
    loop {
        let x = rng.sample(range);
        let y = rng.sample(range);
        let z = rng.sample(range);
        let r2 = x*x + y*y + z*z;
        if r2 <= radius2 {
            return [x, y, z];
        }
    }
}

/// Generating uniformly distributed point on a sphere
/// Alternative method 1 as described in:
/// http://corysimon.github.io/articles/uniformdistn-on-sphere/
pub fn rand_point_on_sphere(radius: f64) -> Point3D {
    debug_assert!(radius > 0.0, "sphere radius cannot be negative: {:?}", radius);

    let mut rng = rand::thread_rng();

    let radius2 = radius*radius;
    let normal = Normal::new(0.0, 10.0);
    // avoid floating point precision lost when dividing by a very small number.
    let min_radius2 = 0.1;
    loop {
        let x = rng.sample(normal);
        let y = rng.sample(normal);
        let z = rng.sample(normal);

        let r2: f64 = x*x + y*y + z*z;
        if r2 > min_radius2 {
            let s = radius/r2.sqrt();
            return [x*s, y*s, z*s];
        }
    }
}

/// Create random points within a sphere
///
/// Parameters
/// ----------
/// * radius: cutoff radius
/// * ntps: the number of points to be generated
pub fn rand_points_within_sphere(radius: f64, npts: usize) -> Points {
    let mut rng = rand::thread_rng();
    let range = Range::new(-radius, radius);

    // using the discarding method, which is simple and also fast
    let radius2 = radius*radius;

    let mut points = vec![];
    loop {
        let x = rng.sample(range);
        let y = rng.sample(range);
        let z = rng.sample(range);
        let r2 = x*x + y*y + z*z;
        if r2 <= radius2 {
            points.push([x, y, z]);
        }
        if points.len() >= npts {
            return points
        }
    }
}
// 9413e1bc-8f8f-4b07-b305-9d2911afabc6 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::976ea5bc-07b5-40f3-bfc8-42e2980d6f31][976ea5bc-07b5-40f3-bfc8-42e2980d6f31]]
pub fn rand_rotate(points: &Points) -> Points {
    let p = rand_point_on_sphere(1.0);
    let v = Vector3::from(p);
    let angle = v.angle(&Vector3::x_axis());
    rotate_about_x_axis(points, angle, [0.0, 0.0, 0.0])
}
// 976ea5bc-07b5-40f3-bfc8-42e2980d6f31 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::4c7139a2-6745-461f-ad04-ea4355283d60][4c7139a2-6745-461f-ad04-ea4355283d60]]
#[test]
fn test_rand_rotate() {
    use std::f64;

    let points = [[-0.02264019, -0.01300713, -0.06295011],
                  [ 1.37326881, -0.01300713, -0.06295011],
                  [-0.44222819, -0.73391213,  0.82834789],
                  [-0.79257355, -1.33584955, -1.69845937],
                  [-0.76587962,  1.29543401, -0.06295011],
                  [-1.46366314,  1.28242565, -0.77914068],
                  [-1.20324889,  1.43034987,  0.82615384]].to_vec();

    let p = rotate_about_x_axis(&points, f64::to_radians(10.0), [0.1, 1.0, 0.9]);
    let expected = [[-0.02264019,  0.16959726, -0.22422758],
                    [ 1.37326881,  0.16959726, -0.22422758],
                    [-0.44222819, -0.69512785,  0.52834576],
                    [-0.79257355, -0.84914501, -2.06459895],
                    [-0.76587962,  1.45816024,  0.00298084],
                    [-1.46366314,  1.56971469, -0.70458806],
                    [-1.20324889,  1.43663514,  0.9020052 ]];

    assert_relative_eq!(p[2][0], expected[2][0], epsilon=1e-4);
    assert_relative_eq!(p[3][1], expected[3][1], epsilon=1e-4);
    assert_relative_eq!(p[4][2], expected[4][2], epsilon=1e-4);

    let rpoints = rand_rotate(&points);

    let npoints = points.len();
    assert_eq!(npoints, rpoints.len());

    for i in 0..npoints {
        for j in (i+1)..npoints {
            let d1 = euclidean_distance(points[i], points[j]);
            let d2 = euclidean_distance(rpoints[i], rpoints[j]);
            assert_relative_eq!(d1, d2, epsilon=1e-4);
        }
    }
}
// 4c7139a2-6745-461f-ad04-ea4355283d60 ends here
