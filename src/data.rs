// [[file:~/Workspace/Programming/gchemol/gchemol.note::fea6623c-f2ad-4d9a-b5d4-8a7c01f7cf01][fea6623c-f2ad-4d9a-b5d4-8a7c01f7cf01]]
use {
    Atom,
    Bond,
    bond::BondKind,
    geometry::euclidean_distance,
};

const COVALENT_RADII: [f64; 96] = [
    0.31,
    0.28,
    1.28,
    0.96,
    0.84,
    0.73,
    0.71,
    0.66,
    0.57,
    0.58,
    1.66,
    1.41,
    1.21,
    1.11,
    1.07,
    1.05,
    1.02,
    1.06,
    2.03,
    1.76,
    1.7,
    1.6,
    1.53,
    1.39,
    1.39,
    1.32,
    1.26,
    1.24,
    1.32,
    1.22,
    1.22,
    1.2,
    1.19,
    1.2,
    1.2,
    1.16,
    2.2,
    1.95,
    1.9,
    1.75,
    1.64,
    1.54,
    1.47,
    1.46,
    1.42,
    1.39,
    1.45,
    1.44,
    1.42,
    1.39,
    1.39,
    1.38,
    1.39,
    1.4,
    2.44,
    2.15,
    2.07,
    2.04,
    2.03,
    2.01,
    1.99,
    1.98,
    1.98,
    1.96,
    1.94,
    1.92,
    1.92,
    1.89,
    1.9,
    1.87,
    1.87,
    1.75,
    1.7,
    1.62,
    1.51,
    1.44,
    1.41,
    1.36,
    1.36,
    1.32,
    1.45,
    1.46,
    1.48,
    1.4,
    1.5,
    1.5,
    2.6,
    2.21,
    2.15,
    2.06,
    2.0,
    1.96,
    1.9,
    1.87,
    1.8,
    1.69
];

pub fn guess_bond_kind(atom1: &Atom, atom2: &Atom) -> BondKind {
    let i1 = atom1.number() - 1;
    let i2 = atom2.number() - 1;
    let cr1 = COVALENT_RADII[i1];
    let cr2 = COVALENT_RADII[i2];
    let d12 = euclidean_distance(atom1.position, atom2.position);

    let rcutoff = (cr1 + cr2)*1.15;
    // println!("{:?}", (i1, i2, rcutoff, d12));
    if d12 > rcutoff {
        return BondKind::Dummy;
    } else {
        return BondKind::Single;
    }
}
// fea6623c-f2ad-4d9a-b5d4-8a7c01f7cf01 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::40a03fbc-5f09-4432-a569-f2216d508de4][40a03fbc-5f09-4432-a569-f2216d508de4]]
include!("data/vdw.txt");
// 40a03fbc-5f09-4432-a569-f2216d508de4 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::0c32acf0-9391-40d5-9c80-3b4bc74f2020][0c32acf0-9391-40d5-9c80-3b4bc74f2020]]
impl Atom {
    /// Access covalent atomic radii
    /// Return None if no data available
    pub fn covalent_radius(&self) -> Option<f64> {
        let i = self.number() - 1;
        if i >= 0 && i< COVALENT_RADII.len() {
            let r = COVALENT_RADII[i];
            Some(r)
        } else {
            None
        }
    }

    pub fn vdw_radius(&self) -> Option<f64> {
        let i = self.number() - 1;
        if i >= 0 && i< VDW_RADII.len() {
            let r = VDW_RADII[i];
            Some(r)
        } else {
            None
        }
    }
}
// 0c32acf0-9391-40d5-9c80-3b4bc74f2020 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::2d6a4ce4-f616-4bea-a5fd-65ed9117e613][2d6a4ce4-f616-4bea-a5fd-65ed9117e613]]
#[test]
fn test_guess_bond_kind() {
    // CH4 molecule
    let atom1 = Atom::new("C", [-0.90203687,  0.62555259,  0.0081889 ]);
    let atom2 = Atom::new("H", [-0.54538244, -0.38325741,  0.0081889 ]);
    let atom3 = Atom::new("H", [-0.54536403,  1.12995078,  0.88184041]);
    let atom4 = Atom::new("H", [-0.54536403,  1.12995078, -0.8654626 ]);
    let atom5 = Atom::new("H", [-1.97203687,  0.62556577,  0.0081889 ]);

    let x = guess_bond_kind(&atom1, &atom2);
    assert!(x != BondKind::Dummy);
    let x = guess_bond_kind(&atom2, &atom3);
    assert!(x == BondKind::Dummy);

    //covalent radii
    assert!(atom1.covalent_radius().is_some());
}
// 2d6a4ce4-f616-4bea-a5fd-65ed9117e613 ends here
