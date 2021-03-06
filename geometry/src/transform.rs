// src

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/geometry/geometry.note::*src][src:1]]
type Point3D = [f64; 3];
type Points = Vec<Point3D>;

use super::prelude::euclidean_distance;

pub use vecfx::nalgebra::{
    Vector3,
    Rotation3,
};

/// Translate all points to a new location
pub fn translate(points: &mut Points, loc: Point3D) {
    for i in 0..points.len() {
        for v in 0..3 {
            points[i][v] += loc[v];
        }
    }
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
// src:1 ends here

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/geometry/geometry.note::*src][src:2]]
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
// src:2 ends here
