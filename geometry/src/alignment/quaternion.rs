// src
// # References
// - B. K. P. Horn, J. Opt. Soc. Am. A, JOSAA, 1987, DOI:10.1364/JOSAA.4.000629.
// - D. L. Theobald, Acta Crystallogr. A, 2005, DOI:10.1107/S0108767305015266.
// - P. Liu, D. K. Agrafiotis, D. L. Theobald, J. Comput. Chem., 2010, DOI:10.1002/jcc.21439.


use crate::base::*;
use nalgebra as na;

pub fn calc_rmsd_rotational_matrix(
    positions_ref: &[[f64; 3]],
    positions_can: &[[f64; 3]],
    weights: Option<&[f64]>,
) -> (f64, [f64; 3], Option<[f64; 9]>) {
    let npts = positions_ref.len();

    // set up weights for atoms
    let default_weights = vec![1.0; npts];
    let weights = weights.unwrap_or(&default_weights);
    // FIXME: Option
    let com_ref = positions_ref.center_of_mass(weights).unwrap();
    let com_can = positions_can.center_of_mass(weights).unwrap();

    // 1. center coordinates of the reference and the candidate
    let mut vectors_ref = positions_ref.to_vectors();
    let mut vectors_can = positions_can.to_vectors();
    for mut v in vectors_ref.iter_mut() {
        *v -= com_ref;
    }
    for mut v in vectors_can.iter_mut() {
        *v -= com_can;
    }

    // 2. computation of the F matrix
    let mut mat_f = na::Matrix3::zeros();
    for i in 0..npts {
        let wi = weights[i];
        mat_f += wi * &vectors_can[i] * &vectors_ref[i].transpose();
    }
    let sxx = mat_f[(0, 0)];
    let sxy = mat_f[(0, 1)];
    let sxz = mat_f[(0, 2)];
    let syx = mat_f[(1, 0)];
    let syy = mat_f[(1, 1)];
    let syz = mat_f[(1, 2)];
    let szx = mat_f[(2, 0)];
    let szy = mat_f[(2, 1)];
    let szz = mat_f[(2, 2)];

    // 3. construct the key matrix K
    let mat_k = na::Matrix4::from_column_slice(&[
        sxx + syy + szz,
        syz - szy,
        szx - sxz,
        sxy - syx,
        syz - szy,
        sxx - syy - szz,
        sxy + syx,
        sxz + szx,
        szx - sxz,
        sxy + syx,
        -sxx + syy - szz,
        syz + szy,
        sxy - syx,
        sxz + szx,
        syz + szy,
        -sxx - syy + szz,
    ]);

    // 4. compute the rotation quaternion
    // which is the eigenvector corresponding to the most positive eigenvalue
    let se = mat_k.symmetric_eigen();
    let ci = se.eigenvalues.imax();
    let q = se.eigenvectors.column(ci);

    // // 5. construct rotation matrix from the quaternion q
    let rot = [
        q[0].powi(2) + q[1].powi(2) - q[2].powi(2) - q[3].powi(2),
        2.0 * (q[1] * q[2] - q[0] * q[3]),
        2.0 * (q[1] * q[3] + q[0] * q[2]),
        2.0 * (q[1] * q[2] + q[0] * q[3]),
        q[0].powi(2) - q[1].powi(2) + q[2].powi(2) - q[3].powi(2),
        2.0 * (q[2] * q[3] - q[0] * q[1]),
        2.0 * (q[1] * q[3] - q[0] * q[2]),
        2.0 * (q[2] * q[3] + q[0] * q[1]),
        q[0].powi(2) - q[1].powi(2) - q[2].powi(2) + q[3].powi(2),
    ];

    // or using nalgebra's library call
    // let q = na::geometry::Quaternion::new(q[0], q[1], q[2], q[3]);
    // let rot = na::geometry::UnitQuaternion::from_quaternion(q).to_rotation_matrix();
    let rotation = Some(rot);

    // 6. calculate superposition rmsd
    let mut rmsd = 0.0f64;
    // rmsd += G_A + G_B
    for i in 0..npts {
        let wi = weights[i];
        let vcan = vectors_can[i] * wi;
        let vref = vectors_ref[i] * wi;
        rmsd += vcan.norm_squared() + vref.norm_squared();
    }

    let emax = se.eigenvalues.as_slice().max();
    let wsum = weights.sum() as f64;
    let rmsd = ((rmsd - 2.0 * emax) / wsum).sqrt();

    // 7. calculate translation
    let mut rotc = [0.0; 3];
    let rotc = if let Some(r) = rotation {
        vec![
            r[0] * com_can[0] + r[1] * com_can[1] + r[2] * com_can[2],
            r[3] * com_can[0] + r[4] * com_can[1] + r[5] * com_can[2],
            r[6] * com_can[0] + r[7] * com_can[1] + r[8] * com_can[2],
        ]
    } else {
        com_can.as_slice().to_vec()
    };

    let trans = [
        rotc[0] - com_ref[0],
        rotc[1] - com_ref[1],
        //rotc[2] - com_ref[2],
        dbg!(rotc[2]) - dbg!(com_ref[2]),
    ];

    return (rmsd, trans, rotation);
}

// test
// This is an example where the original QCP program gives the incorrect result.

#[test]
fn test_quaterion() {
    use super::qcprot;
    use gchemol::Molecule;
    use gchemol::prelude::*;

    let positions_ref = vec![[ 0.     ,  0.     ,  0.50222],
                             [ 0.     ,  0.     ,  1.57133],
                             [ 0.     ,  0.     , -0.65495]];

    let positions_can = vec![[-0.35846,  0.02588,  0.47287],
                             [ 0.28996,  1.02611,  0.47287],
                             [ 0.0685 , -1.05199,  0.47287]];

    let (r1, tran1, rot1) = calc_rmsd_rotational_matrix(&positions_ref, &positions_can, None);
    let (r2, tran2, rot2) = qcprot::calc_rmsd_rotational_matrix(&positions_ref, &positions_can, None);
    assert_relative_eq!(r1, r2, epsilon=1e-3);
    dbg!(tran1);
    dbg!(tran2);

    for i in 0..3 {
        assert_relative_eq!(tran1[i], tran2[i], epsilon=1e-3);
    }
    let rot1 = rot1.expect("rot1");
    let rot2 = rot2.expect("rot2");
    for i in 0..9 {
        assert_relative_eq!(rot1[i], rot2[i], epsilon=1e-3);
    }
}
