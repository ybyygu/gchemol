// basic stats

#[macro_use]
extern crate approx;

use gchemol::prelude::*;

#[test]
fn test_stats() {
    let p = [-1.0, 1.0, 5.2];

    // provides common statistics methods
    // for more: https://doc.rust-lang.org/test/stats/trait.Stats.html
    assert_eq!(5.2, p.max());
    assert_eq!(-1.0, p.min());
    assert_relative_eq!(1.7333, p.mean(), epsilon=1e-4);
}

// 3D points

#[test]
fn test_point3() {
    // methods for points in 3D space
    let p1 = [ 0.471321,  0.430787,  0.      ];
    let p2 = [-0.218686, -0.495915,  0.      ];

    assert_relative_eq!(1.1554, p1.distance(p2), epsilon=1e-4);
}

// 3D positions in cartesian coordinates

#[test]
fn test_positions() {
    // construct array
    let positions = [[ 0.471321,  0.430787,  0.      ],
                     [-0.253529,  1.463036,  0.      ],
                     [-0.218686, -0.495915,  0.      ]];

    // as flat 1D slice
    let positions_flat = positions.as_flat();
    assert_eq!(positions_flat.len(), 9);
    assert_eq!(-0.253529, positions_flat[3]);

    // norm
    assert_relative_eq!(1.7048, positions.norm(), epsilon=1e-4);
    assert_relative_eq!(positions_flat.norm(), positions.norm(), epsilon=1e-4);

    // convert to nalgebra 3xN matrix (column major)
    let mat = positions.to_dmatrix();

    // matrix dot
    let mat2 = &mat * &mat;
    let expected = [[ 0.11292649,  0.83329585,  0.        ],
                    [-0.4904156 ,  2.03125734,  0.        ],
                    [ 0.02265753, -0.81974858,  0.        ]].to_dmatrix();

    assert_relative_eq!(mat2, expected, epsilon=1e-4);

    // as 1D slice
    let slice = mat2.as_slice();
    assert_eq!(9, slice.len());

    // as 3x3 array
    let arr = slice.as_positions();
    assert_eq!(arr.len(), 3);

    // mutable version
    // construct vector
    let mut positions = vec![[ 0.471321,  0.430787,  0.      ],
                             [-0.253529,  1.463036,  0.      ],
                             [-0.218686, -0.495915,  0.      ]];

    let positions_flat = positions.as_mut_flat();
    positions_flat[3] = 0.0;
    assert_eq!(0.0, positions_flat[3]);
}

// TODO alignment


