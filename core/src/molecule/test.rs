// [[file:~/Workspace/Programming/gchemol/core/gchemol-core.note::*test.rs][test.rs:1]]
use super::*;

#[test]
fn test_molecule_basic() {
    // construct molecule
    let mut mol = Molecule::new("test");
    assert_eq!("test", mol.name);

    let mut mol = Molecule::default();
    let atom1 = Atom::new("Fe", [1.2; 3]);
    let atom2 = Atom::new("Fe", [1.0; 3]);
    let atom3 = Atom::new("C", [0.0; 3]);
    let atom4 = Atom::new("O", [2.1; 3]);
    let a1 = mol.add_atom(atom1);
    let a2 = mol.add_atom(atom2);
    let a3 = mol.add_atom(atom3);
    let a4 = mol.add_atom(atom4);
    assert_eq!(4, mol.natoms());

    let b1 = mol.add_bond(a1, a2, Bond::default());
    let b2 = mol.add_bond(a3, a4, Bond::default());
    assert_eq!(2, mol.nbonds());
    mol.remove_bond_between(a1, a2).expect("failed to remove bond between a1 and a2");
    mol.remove_bond(b2).expect("failed to remove bond b2");
    assert_eq!(0, mol.nbonds());
    let b14 = mol.add_bond(a1, a4, Bond::default());
    assert_eq!(1, mol.nbonds());

    // bonded partners
    let real_b14 = mol.get_bond(b14).expect("failed to get bond b14");
    assert_eq!(real_b14.index(), b14);
    let (n1, n4) = mol.partners(&b14).expect("failed to get bond partners using bond index");
    assert_eq!(n1.index(), a1.index());
    assert_eq!(n4.index(), a4.index());
    // get partners using bond
    let (n1, n4) = real_b14.partners(&mol).expect("failed to get bond partners using bond struct");
    assert_eq!(n1.index(), a1);
    assert_eq!(n4.index(), a4);

    // get atom neighbors
    let indices = mol.neighbors(a1);
    assert_eq!(1, indices.len());
    assert!(indices.contains(&a4));

    // loop over atoms
    for a in mol.atoms() {
        //
    }

    // loop over bonds
    for b in mol.bonds() {
        //
    }

    // pick a single atom
    let a = mol.get_atom(AtomIndex::new(0)).expect("failed to get atom with index 0");
    assert_eq!("Fe", a.symbol());
    assert_eq!(1.2, a.position()[0]);
    let a = mol.get_atom(a1).expect("failed to get atom a1");
    assert_eq!("Fe", a.symbol());
    assert_eq!(1.2, a.position()[0]);
}

#[test]
fn test_molecule_scaled_position() {
    let mut mol = Molecule::default();
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    let lat = Lattice::default();
    mol.set_lattice(lat);
    let fxyzs = [[ 0.41737596,  0.75597855,  0.16098257],
                 [ 0.41565679,  0.68546917,  0.19264617],
                 [ 0.38461882,  0.68391421,  0.22795131],
                 [ 0.37458942,  0.74128686,  0.24842385]].to_vec();

    mol.set_scaled_positions(&fxyzs).unwrap();
    let fxyzs2 = mol.scaled_positions().unwrap();
    assert_eq!(fxyzs, fxyzs2);
}

#[test]
fn test_molecule_other() {
    let mut mol = Molecule::default();
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    mol.add_atom(Atom::default());
    //set atom positions
    let positions = [[-0.90203687,  0.62555259,  0.0081889 ],
                     [-0.54538244, -0.38325741,  0.0081889 ],
                     [-0.54536403,  1.12995078, -0.8654626 ],
                     [-1.97203687,  0.62556577,  0.0081889 ]];
    mol.set_positions(&positions).unwrap();
    let a = mol.get_atom(AtomIndex::new(0)).expect("failed to get atom with index 0");
    assert_eq!(a.position()[0], -0.90203687);

    // loop over fragments
    let frags = mol.fragment();
    assert_eq!(mol.nfragments(), frags.len());
    for m in frags {
        m.formula();
    }
}

#[test]
fn test_molecule_rebond() {
    let atom1 = Atom::new("C", [-0.90203687,  0.62555259,  0.0081889 ]);
    let atom2 = Atom::new("H", [-0.54538244, -0.38325741,  0.0081889 ]);
    let atom3 = Atom::new("H", [-0.54536403,  1.12995078,  0.88184041]);
    let atom4 = Atom::new("H", [-0.54536403,  1.12995078, -0.8654626 ]);
    let atom5 = Atom::new("H", [-1.97203687,  0.62556577,  0.0081889 ]);

    let mut mol = Molecule::default();
    mol.add_atom(atom1);
    mol.add_atom(atom2);
    mol.add_atom(atom3);
    mol.add_atom(atom4);
    mol.add_atom(atom5);

    assert_eq!(5, mol.natoms());
    assert_eq!(0, mol.nbonds());
    mol.rebond();
    assert_eq!(4, mol.nbonds());
}
// test.rs:1 ends here
