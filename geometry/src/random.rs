// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::9413e1bc-8f8f-4b07-b305-9d2911afabc6][9413e1bc-8f8f-4b07-b305-9d2911afabc6]]
use rand::{self, Rng};
use rand::distributions::{Range, Normal};

use super::transform::*;
use super::prelude::euclidean_distance;

type Point3D = [f64; 3];
type Points = Vec<Point3D>;

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

// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::976ea5bc-07b5-40f3-bfc8-42e2980d6f31][976ea5bc-07b5-40f3-bfc8-42e2980d6f31]]
pub fn rand_rotate(points: &Points) -> Points {
    let p = rand_point_on_sphere(1.0);
    let v = Vector3::from(p);
    let angle = v.angle(&Vector3::x_axis());
    rotate_about_x_axis(points, angle, [0.0, 0.0, 0.0])
}
// 976ea5bc-07b5-40f3-bfc8-42e2980d6f31 ends here

// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::4c7139a2-6745-461f-ad04-ea4355283d60][4c7139a2-6745-461f-ad04-ea4355283d60]]
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
