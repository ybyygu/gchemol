// [[file:~/Workspace/Programming/gchemol/gchemol.note::6746df7d-0871-43bf-98cd-2e00e15020a5][6746df7d-0871-43bf-98cd-2e00e15020a5]]
use cgmath::{Vector3, Matrix3, Point3, Deg};
use cgmath::prelude::*;

fn cart_to_frac(matrix: Matrix3<f64>,
                coordinates: Vec<Vector3<f64>>) -> Vec<Vector3<f64>>
{
    let mut fractional = Vec::new();
    let inv = matrix.transpose().invert().unwrap();
    for v in coordinates {
        fractional.push(inv*v);
    }

    fractional
}
// 6746df7d-0871-43bf-98cd-2e00e15020a5 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c57f4ca0-4e68-4e30-a91b-7cbd47b7071c][c57f4ca0-4e68-4e30-a91b-7cbd47b7071c]]
use std::f64;

fn get_nearest_image(
    cell: Matrix3<f64>,
    position1: Point3<f64>,
    position2: Point3<f64>) -> (Vector3<f64>, f64)
{
    let d = position1.distance(position2);

    // loop 27 possible point images
    let relevant_images = [-1, 0, 1];
    let mut distance = f64::MAX;
    let mut image = Vector3::from_value(0_f64);
    for x in relevant_images.iter() {
        for y in relevant_images.iter() {
            for z in relevant_images.iter() {
                let p = position2 + (*x as f64)*cell.x + (*y as f64)*cell.y + (*z as f64)*cell.z;
                let d = position1.distance(p);
                if d < distance {
                    distance = d;
                    image.x = *x as f64;
                    image.y = *y as f64;
                    image.z = *z as f64;
                }
            }
        }
    }

    (image, distance)
}

#[test]
fn test_get_nearest_image() {
    let mat1 = Matrix3::new(5.09, 0.00, 0.00,
                            0.00, 6.74, 0.00,
                            0.00, 0.00, 4.53);

    let p1  = Point3::new(0.18324000,   1.68500000,   3.85050000);
    let p13 = Point3::new(4.53010000,   1.68500000,   2.03850000);
    let p10 = Point3::new(0.94674000,   2.94538000,   1.48584000);
    let dp1_13 = 1.95847;
    let dp1_10 = 2.61920;

    let (image, d) = get_nearest_image(mat1, p1, p13);
    assert_relative_eq!(d, dp1_13, epsilon=1e-4);
    assert_relative_eq!(image.x, -1.0, epsilon=1e-4);
    assert_relative_eq!(image.y, 0.0, epsilon=1e-4);
    assert_relative_eq!(image.z, 0.0, epsilon=1e-4);

    let (image, d) = get_nearest_image(mat1, p1, p10);
    assert_relative_eq!(d, dp1_10, epsilon=1e-4);
}
// c57f4ca0-4e68-4e30-a91b-7cbd47b7071c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::c90b9d01-096f-47c1-bbfd-2649862e61dc][c90b9d01-096f-47c1-bbfd-2649862e61dc]]
fn cell_vectors_to_parameters(matrix: Matrix3<f64>) -> (f64, f64, f64, f64, f64, f64) {
    let a = matrix.x.magnitude();
    let b = matrix.y.magnitude();
    let c = matrix.z.magnitude();

    let alpha: Deg<_> = matrix.y.angle(matrix.z).into();
    let beta: Deg<_> = matrix.x.angle(matrix.z).into();
    let gamma: Deg<_> = matrix.x.angle(matrix.y).into();

    (a, b, c, alpha.0, beta.0, gamma.0)
}
// c90b9d01-096f-47c1-bbfd-2649862e61dc ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::9cab3b07-9781-48cd-a6fb-e6ee248b93dd][9cab3b07-9781-48cd-a6fb-e6ee248b93dd]]
#[test]
fn test_cell() {
    // ovito/tests/files/LAMMPS/multi_sequence_1.dump
    let mat1 = Matrix3::new(5.09, 0.00, 0.00,
                            0.00, 6.74, 0.00,
                            0.00, 0.00, 4.53);
    let inv = mat1.transpose().invert().unwrap();

    let v1 = Vector3::new(2.1832, 1.6850, 3.8505);
    let v2 = Vector3::new(6.9068, 5.0550, 0.6795);
    let v3 = Vector3::new(4.3618, 5.0550, 1.5855);

    let fracs = cart_to_frac(mat1, vec![v1, v2, v3]);
    assert_relative_eq!(fracs[0].x, 0.4289, epsilon=1e-3);
    assert_relative_eq!(fracs[0].y, 0.2500, epsilon=1e-3);
    assert_relative_eq!(fracs[0].z, 0.8500, epsilon=1e-3);
    assert_relative_eq!(fracs[1].x, 1.3569, epsilon=1e-3);
    assert_relative_eq!(fracs[2].z, 0.3500, epsilon=1e-3);

    let mat2 = Matrix3::new(15.3643, 0.0, 0.0,
                            4.5807, 15.5026, 0.0,
                            0.0, 0.0, 17.4858);

    let (a, b, c, alpha, beta, gamma) = cell_vectors_to_parameters(mat2);
    assert_relative_eq!(a, 15.3643, epsilon=1e-4);
    assert_relative_eq!(b, 16.1652, epsilon=1e-4);
    assert_relative_eq!(c, 17.4858, epsilon=1e-4);

    assert_relative_eq!(alpha, 90.0, epsilon=1e-4);
    assert_relative_eq!(beta, 90.0, epsilon=1e-4);
    assert_relative_eq!(gamma, 73.5386, epsilon=1e-4);
}
// 9cab3b07-9781-48cd-a6fb-e6ee248b93dd ends here
