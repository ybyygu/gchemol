// [[file:~/Workspace/Programming/gchemol/gchemol.note::53cbd3c0-e164-4bad-b535-6fd6df916650][53cbd3c0-e164-4bad-b535-6fd6df916650]]
#![feature(conservative_impl_trait)]
#![allow(dead_code)]

#[macro_use]
extern crate error_chain;

// We'll put our errors in an `errors` module, and other modules in
// this crate will `use errors::*;` to get access to everything
// `error_chain!` creates.
pub mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! { }
}

extern crate petgraph;
extern crate cgmath;
#[macro_use]
extern crate timeit;
#[macro_use]
extern crate approx;
extern crate rand;
extern crate itertools;

pub type Point3D = [f64; 3];
pub type Points = Vec<Point3D>;

pub mod element;
pub use element::AtomKind;
pub use element::AtomKind::{Element, Dummy};
pub mod atom;
pub use atom::Atom;
pub mod bond;
pub use bond::Bond;
pub mod graph;
pub mod molecule;
pub mod topology;
pub mod io;
pub use io::write_as_xyz;

#[inline]
pub fn euclidean_distance(p1: Point3D, p2: Point3D) -> f64 {
    let mut d2 = 0.0;
    for v in 0..3 {
        let dv = p2[v] - p1[v];
        d2 += dv*dv;
    }

    d2.sqrt()
}
// 53cbd3c0-e164-4bad-b535-6fd6df916650 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::0cd950e1-e3e1-44cf-96f8-8629f3bf3f95][0cd950e1-e3e1-44cf-96f8-8629f3bf3f95]]
use rand::Rng;
use rand::distributions::{Distribution, Range, Normal};

use itertools::Itertools;

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
// 0cd950e1-e3e1-44cf-96f8-8629f3bf3f95 ends here

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
        for j in i..npoints {
            let d1 = euclidean_distance(points[i], points[j]);
            let d2 = euclidean_distance(rpoints[i], rpoints[j]);
            assert_relative_eq!(d1, d2, epsilon=1e-4);
        }
    }
}
// 46344373-08ce-4242-bf42-5981f2ee1bd0 ends here
