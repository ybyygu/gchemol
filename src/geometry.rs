// [[file:~/Workspace/Programming/gchemol/gchemol.note::f2ed012d-0e00-4288-b59e-0cb61f7921c2][f2ed012d-0e00-4288-b59e-0cb61f7921c2]]
use Point3D;
use Points;

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

// [[file:~/Workspace/Programming/gchemol/gchemol.note::9413e1bc-8f8f-4b07-b305-9d2911afabc6][9413e1bc-8f8f-4b07-b305-9d2911afabc6]]
use rand::{self, Rng};
use rand::distributions::{Distribution, Range, Normal};

/// create a random point within a sphere
/// References
/// ----------
/// https://stackoverflow.com/a/5408344
pub fn rand_point_within_sphere(radius: f64) -> Point3D {
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
// 9413e1bc-8f8f-4b07-b305-9d2911afabc6 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::46344373-08ce-4242-bf42-5981f2ee1bd0][46344373-08ce-4242-bf42-5981f2ee1bd0]]
use cgmath::prelude::*;
use cgmath::{Quaternion, Vector3};

/// generate a random quaternion for rotation
fn rand_quaternion() -> Quaternion<f64> {
    let radius = 1.0;
    let p = rand_point_on_sphere(radius);
    let v = Vector3::from(p);
    let s = (v.magnitude2() + radius*radius).sqrt();

    Quaternion::from_sv(radius/s, v/s)
}

pub fn rand_rotate(points: Points) -> Points {
    let rot = rand_quaternion();
    let mut rpoints = vec![];
    for &p in points.iter() {
        let v = Vector3::from(p);
        let t: Point3D = (rot*v).into();
        rpoints.push(t);
    }

    rpoints
}

#[test]
fn test_rand_rotate() {
    let points = [[-0.02264019, -0.01300713, -0.06295011],
                  [ 1.37326881, -0.01300713, -0.06295011],
                  [-0.44222819, -0.73391213,  0.82834789],
                  [-0.79257355, -1.33584955, -1.69845937],
                  [-0.76587962,  1.29543401, -0.06295011],
                  [-1.46366314,  1.28242565, -0.77914068],
                  [-1.20324889,  1.43034987,  0.82615384]];

    let rpoints = rand_rotate(points.to_vec());

    let npoints = points.len();
    assert_eq!(npoints, rpoints.len());
    assert!(points[0][0] != rpoints[0][0]);

    for i in 0..npoints {
        for j in (i+1)..npoints {
            let d1 = euclidean_distance(points[i], points[j]);
            let d2 = euclidean_distance(rpoints[i], rpoints[j]);
            assert_relative_eq!(d1, d2, epsilon=1e-4);
        }
    }
}
// 46344373-08ce-4242-bf42-5981f2ee1bd0 ends here
