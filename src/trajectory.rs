// [[file:~/Workspace/Programming/gchemol/gchemol.note::e2130a32-e39f-4b7b-9014-515f18ff5f48][e2130a32-e39f-4b7b-9014-515f18ff5f48]]
use molecule::Molecule;

impl Molecule {
    /// Interpolate between this structure and end_structure. Useful for
    /// construction of NEB inputs.
    fn interpolate(&self, end_mol: &Molecule, nimages: usize) -> Vec<Molecule> {
        unimplemented!()
    }
}
// e2130a32-e39f-4b7b-9014-515f18ff5f48 ends here
