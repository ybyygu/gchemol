// io.rs
// :PROPERTIES:
// :header-args: :tangle tests/io.rs :comments org
// :ID:       75afbab7-abee-46f6-ae31-1232b81842ce
// :END:
// demonstrate how to read/write molecules

extern crate gchemol;
extern crate tempfile;

#[test]
fn io_from_file_to_file() {
    use gchemol::Molecule;
    use gchemol::prelude::*;

    // 1. string io
    let txt1 = String::from_file("tests/files/mol2/alanine-gv.mol2").expect("txt1 from file");

    // write to a temp file
    let tfile = tempfile::Builder::new()
        .suffix(".mol2")
        .tempfile().expect(".mol2 temp file.");

    // validate
    txt1.to_file(&tfile).expect("write mol2");
    let txt2 = String::from_file(tfile).expect("txt2 from file");
    assert_eq!(txt1, txt2);

    // 2. molecule io
    let mol1 = Molecule::from_file("tests/files/mol2/alanine-gv.mol2").expect("mol1 from file");

    // write to a temp file in xyz format
    let tfile = tempfile::Builder::new()
        .suffix(".xyz")
        .tempfile().expect(".mol2 temp file.");

    mol1.to_file(&tfile).expect("write mol2 file");
    let mol2 = Molecule::from_file(tfile).expect("mol2 from file");
    assert_eq!(mol1.natoms(), mol2.natoms());
}
