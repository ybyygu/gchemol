// [[file:~/Workspace/Programming/gchemol/geometry/geometry.note::8a489c41-809c-46c9-8a23-5d926db1c5ca][8a489c41-809c-46c9-8a23-5d926db1c5ca]]
use super::base::*;
use quicli::prelude::*;

// direct linear interpolation
// # Parameters
// - positions_initial: initial coordinates
// - positions_final  : final coordinates
// - nimages (> 2): total number of images to be generated (including two endpoints)
//
fn linear_interpolate_positions(
    positions_initial : &[[f64; 3]],
    positions_final   : &[[f64; 3]],
    nimages           : usize) -> Result<Vec<Vec<[f64; 3]>>>
{
    if nimages <= 2 {
        bail!("The number of images is less than 2 (endpoints included).");
    }

    let npts = positions_initial.len();
    if npts != positions_final.len() {
        bail!("Different number of points!");
    }

    let im = positions_initial.to_dmatrix();
    let fm = positions_final.to_dmatrix();
    let n = nimages - 1;
    let dv = (&fm - &im) / n as f64;

    let mut images = Vec::with_capacity(nimages);
    for i in 0..nimages {
        let p = &im + &dv * i as f64;
        let mut image = Vec::with_capacity(npts);
        for j in 0..npts {
            let mut cj = p.column(j).clone_owned();
            let xyz: [f64; 3] = cj.into();
            image.push(xyz);
        }
        images.push(image);
    }

    Ok(images)
}

#[test]
fn test_interpolate() {
    let ipositions = vec![[-5.102629,  1.482321,  0.      ],
                          [-5.702629,  1.482321,  0.      ]];

    let fpositions = vec![[-4.737629,  1.482321,  0.      ],
                          [-6.067629,  1.482321,  0.      ]];

    let x = linear_interpolate_positions(&ipositions, &fpositions, 3).unwrap();
    assert_eq!(3, x.len());
}
// 8a489c41-809c-46c9-8a23-5d926db1c5ca ends here
