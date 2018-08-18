// [[file:~/Workspace/Programming/gchemol/geometry.note::d5604dc1-f9b1-4dca-a3c0-199cdd4ec28f][d5604dc1-f9b1-4dca-a3c0-199cdd4ec28f]]
use geometry::prelude::*;
use nalgebra as na;

fn quaternion_rotate(positions_ref: &[[f64; 3]], positions_can: &[[f64; 3]]) {
    let npts = positions_ref.len();
    let cog_ref = positions_ref.center_of_geometry();
    let cog_can = positions_can.center_of_geometry();

    // center coordinates
    let mut vectors_ref = positions_ref.to_vectors();
    let mut vectors_can = positions_can.to_vectors();

    for mut v in vectors_ref.iter_mut() {
        *v -= cog_ref;
    }
    for mut v in vectors_can.iter_mut() {
        *v -= cog_can;
    }

    let pm1 = Positions::from_columns(&vectors_ref);
    let pm2 = Positions::from_columns(&vectors_can);

    // Computation of the covariance matrix
    let mat_cov = &pm2 * &pm1.transpose();

    let r11 = mat_cov[(0, 0)];
    let r12 = mat_cov[(0, 1)];
    let r13 = mat_cov[(0, 2)];
    let r21 = mat_cov[(1, 0)];
    let r22 = mat_cov[(1, 1)];
    let r23 = mat_cov[(1, 2)];
    let r31 = mat_cov[(2, 0)];
    let r32 = mat_cov[(2, 1)];
    let r33 = mat_cov[(2, 2)];

    // compute the rotation quaternion
    println!("{}", mat_cov);
    let f = [
        r11 + r22 + r33, r23 - r32, r31 - r13, r12 - r21,
        r23 - r32, r11 - r22 - r33, r12 + r21, r13 + r31,
        r31 - r13, r12 + r21, -r11 + r22 - r33, r23 + r32,
        r12 - r21, r13 + r31, r23 + r32, -r11 - r22 + r33
    ];
    let mat_f = na::Matrix4::<f64>::from_column_slice(&f);

    // eigenvector corresponding to the most positive eigenvalue
    let se = mat_f.symmetric_eigen();
    let ci = se.eigenvalues.imax();
    let vq = se.eigenvectors.column(ci);

    // construct rotation matrix from the quaternion q
    let [q0, q1, q2, q3] = [vq[0], vq[1], vq[2], vq[3]];

    let r_q = [
        q0.powi(2) + q1.powi(2) - q2.powi(2) - q3.powi(2),
        2.0 * (q1 * q2 - q0 * q3),
        2.0 * (q1 * q3 + q0 * q2),
        2.0 * (q1 * q2 + q0 * q3),
        q0.powi(2) - q1.powi(2) + q2.powi(2) - q3.powi(2),
        2.0 * (q2 * q3 - q0 * q1),
        2.0 * (q1 * q3 - q0 * q2),
        2.0 * (q2 * q3 + q0 * q1),
        q0.powi(2) - q1.powi(2) - q2.powi(2) + q3.powi(2)];

    let mat_r = na::Matrix3::<f64>::from_column_slice(&r_q);

    let mut positions_new = vec![];
    for mut v in vectors_can.iter() {
        let p = (v.transpose() * mat_r).transpose() + cog_ref;
        positions_new.push(p);
    }

    println!("{:#?}", positions_new);

    // calculate rmsd
    let emax = se.eigenvalues.max();

    let mut rmsd = 0.0f64;
    for i in 0..npts {
        let vcan = vectors_can[i];
        let vref = vectors_ref[i];
        rmsd += vcan.norm_squared() + vref.norm_squared() - 2.0 * emax;
    }

    rmsd /= npts as f64;
    let rmsd = rmsd.sqrt();
    println!("rmsd = {:6.3}", rmsd);
}
// d5604dc1-f9b1-4dca-a3c0-199cdd4ec28f ends here
