// [[file:~/Workspace/Programming/gchemol/gchemol.note::4b054000-d421-4fb6-a18f-5379838d89fc][4b054000-d421-4fb6-a18f-5379838d89fc]]
use Point3D;
use Points;
use euclidean_distance;

/// check if any pair of points come too close
fn close_contact(points: &Points) -> bool {
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
fn get_distance_matrix(points: Points) -> Vec<Vec<f64>>{
    let npts = points.len();

    // fill distance matrix
    let mut distmat = vec![];
    for i in 0..npts {
        let mut dijs = vec![];
        for j in i..npts {
            let dij = euclidean_distance(points[i], points[j]);
            dijs.push(dij);
        }
        distmat.push(dijs);
    }

    debug_assert!(distmat.len() == npts, "distance matrix shape is incorrect!");
    debug_assert!(distmat[0].len() == npts, "distance matrix shape is incorrect!");
    distmat
}
// 4b054000-d421-4fb6-a18f-5379838d89fc ends here
