// interplate

// [[file:~/Workspace/Programming/gchemol/core/gchemol-core.note::*interplate][interplate:1]]
use crate::molecule::Molecule;

impl Molecule {
    /// Interpolate between this structure and end_structure. Useful for
    /// construction of NEB inputs.
    fn interpolate(&self, end_mol: &Molecule, nimages: usize) -> Vec<Molecule> {
        unimplemented!()
    }
}
// interplate:1 ends here
