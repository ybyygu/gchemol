// [[file:~/Workspace/Programming/gchemol/gchemol.note::7e391e0e-a3e8-4c22-b881-e0425d0926bc][7e391e0e-a3e8-4c22-b881-e0425d0926bc]]
//===============================================================================#
//   DESCRIPTION:  molecule repsented in graph data structure
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  based on petgraph
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-12 Thu 15:48>
//       UPDATED:  <2018-08-06 Mon 11:18>
//===============================================================================#

use std::collections::HashMap;
use petgraph;
use petgraph::prelude::*;
use std::convert;
use quicli::prelude::*;

pub use {
    Point3D,
    Points,
    lattice::Lattice,
};

mod view;
mod formula;
mod fragment;
mod property;
mod clean;
mod rebond;
mod order;
mod edit;

#[cfg(test)]
mod test;

use self::property::*;
// 7e391e0e-a3e8-4c22-b881-e0425d0926bc ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a806642c-37da-4ce1-aa7b-0fb8d00233e3][a806642c-37da-4ce1-aa7b-0fb8d00233e3]]
use std::fmt::{self, Debug, Display};

const ELEMENT_DATA: [(&'static str, &'static str); 118] = [
    ("H", "hydrogen"),
    ("He", "helium"),
    ("Li", "lithium"),
    ("Be", "beryllium"),
    ("B", "boron"),
    ("C", "carbon"),
    ("N", "nitrogen"),
    ("O", "oxygen"),
    ("F", "fluorine"),
    ("Ne", "neon"),
    ("Na", "sodium"),
    ("Mg", "magnesium"),
    ("Al", "aluminum"),
    ("Si", "silicon"),
    ("P", "phosphorus"),
    ("S", "sulfur"),
    ("Cl", "chlorine"),
    ("Ar", "argon"),
    ("K", "potassium"),
    ("Ca", "calcium"),
    ("Sc", "scandium"),
    ("Ti", "titanium"),
    ("V", "vanadium"),
    ("Cr", "chromium"),
    ("Mn", "manganese"),
    ("Fe", "iron"),
    ("Co", "cobalt"),
    ("Ni", "nickel"),
    ("Cu", "copper"),
    ("Zn", "zinc"),
    ("Ga", "gallium"),
    ("Ge", "germanium"),
    ("As", "arsenic"),
    ("Se", "selenium"),
    ("Br", "bromine"),
    ("Kr", "krypton"),
    ("Rb", "rubidium"),
    ("Sr", "strontium"),
    ("Y", "yttrium"),
    ("Zr", "zirconium"),
    ("Nb", "niobium"),
    ("Mo", "molybdenum"),
    ("Tc", "technetium"),
    ("Ru", "ruthenium"),
    ("Rh", "rhodium"),
    ("Pd", "palladium"),
    ("Ag", "silver"),
    ("Cd", "cadmium"),
    ("In", "indium"),
    ("Sn", "tin"),
    ("Sb", "antimony"),
    ("Te", "tellurium"),
    ("I", "iodine"),
    ("Xe", "xenon"),
    ("Cs", "cesium"),
    ("Ba", "barium"),
    ("La", "lanthanum"),
    ("Ce", "cerium"),
    ("Pr", "praesodymium"),
    ("Nd", "neodymium"),
    ("Pm", "promethium"),
    ("Sm", "samarium"),
    ("Eu", "europium"),
    ("Gd", "gadolinium"),
    ("Tb", "terbium"),
    ("Dy", "dyprosium"),
    ("Ho", "holmium"),
    ("Er", "erbium"),
    ("Tm", "thulium"),
    ("Yb", "ytterbium"),
    ("Lu", "lutetium"),
    ("Hf", "hafnium"),
    ("Ta", "tantalium"),
    ("W", "wolfram"),
    ("Re", "rhenium"),
    ("Os", "osmium"),
    ("Ir", "iridium"),
    ("Pt", "platinum"),
    ("Au", "gold"),
    ("Hg", "mercury"),
    ("Tl", "thallium"),
    ("Pb", "lead"),
    ("Bi", "bismuth"),
    ("Po", "polonium"),
    ("At", "astatine"),
    ("Rn", "radon"),
    ("Fr", "francium"),
    ("Ra", "radium"),
    ("Ac", "actinium"),
    ("Th", "thorium"),
    ("Pa", "protactinium"),
    ("U", "uranium"),
    ("Np", "neptunium"),
    ("Pu", "plutonium"),
    ("Am", "americium"),
    ("Cm", "curium"),
    ("Bk", "berkelium"),
    ("Cf", "californium"),
    ("Es", "einsteinium"),
    ("Fm", "fermium"),
    ("Mv", "mendelevium"),
    ("No", "nobelium"),
    ("Lr", "lawrencium"),
    ("Rf", "rutherfordium"),
    ("Db", "dubnium"),
    ("Sg", "seaborgium"),
    ("Bh", "bohrium"),
    ("Hs", "hassium"),
    ("Mt", "meitnerium"),
    ("Uun", "ununnilium"),
    ("Uuu", "unununium"),
    ("Uub", "ununbium"),
    ("Uut", "ununtrium"),
    ("Uuq", "ununquadium"),
    ("Uup", "ununpentium"),
    ("Uuh", "ununhexium"),
    ("Uus", "ununseptium"),
    ("Uuo", "ununoctium")];
// a806642c-37da-4ce1-aa7b-0fb8d00233e3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::731651b9-ebba-4df3-86f7-083f837e4065][731651b9-ebba-4df3-86f7-083f837e4065]]
#[derive(Debug, Clone, PartialEq)]
pub enum AtomKind {
    /// physical elements
    Element(usize),
    /// a dummy-atom is not a real atom
    Dummy(String),
}

impl AtomKind {
    pub fn symbol(&self) -> &str {
        match self {
            &Element(num) => ELEMENT_DATA[num-1].0,
            &Dummy(ref sym) => sym
        }
    }

    pub fn number(&self) -> usize {
        match self {
            &Element(num) => num,
            &Dummy(ref sym) => 0,
        }
    }

    pub fn name(&self) -> String {
        match self {
            &Element(num) => String::from(ELEMENT_DATA[num-1].1),
            &Dummy(ref sym) => format!("dummy atom {}", sym),
        }
    }
}

impl fmt::Display for AtomKind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            &Element(num) => write!(f, "{:}", self.symbol()),
            &Dummy(ref sym) =>  write!(f, "{:}", self.symbol()),
        }
    }
}
// 731651b9-ebba-4df3-86f7-083f837e4065 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b95edc21-e696-4625-ba99-94257394772d][b95edc21-e696-4625-ba99-94257394772d]]
use self::AtomKind::{Element, Dummy};

/// Return AtomKind using common sense
pub fn atom_kind_from_string<T: Into<String>>(sym: T) -> AtomKind {
    let sym = sym.into();

    // element specified in number
    if let Ok(x) = sym.parse::<usize>() {
        return Element(x);
    }

    // element specified in symbol or long name
    for (i, &(s, n) ) in ELEMENT_DATA.iter().enumerate() {
        if s == sym  || n == sym {
            return Element(i+1);
        }
    }

    // treat as dummy atom for the last resort
    Dummy(sym)
}

#[test]
fn test_element() {
    let x = Element(12);
    assert_eq!(12, x.number());
    assert_eq!("Mg", x.symbol());
    assert_eq!("magnesium", x.name());

    let x = Dummy("X".to_string());
    assert_eq!("X", x.symbol());
    assert_eq!(0, x.number());

    let k = atom_kind_from_string("11");
    assert_eq!(k.number(), 11);
    assert_eq!(k.symbol(), "Na");
    assert_eq!(k.name(), "sodium");
}
// b95edc21-e696-4625-ba99-94257394772d ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::150189fd-57d9-4e19-a888-d64497f5ba7e][150189fd-57d9-4e19-a888-d64497f5ba7e]]
use std::hash::Hash;
use std::cmp::Ordering;

#[derive (Debug, Clone)]
/// Atom is the smallest particle still characterizing a chemical element.
///
/// # Reference
///
/// https://goldbook.iupac.org/html/A/A00493.html
///
pub struct Atom {
    /// Arbitrary property stored in key-value pair.
    /// Key is a string type, but it is the responsibility
    /// of the user to correctly interpret the value.
    pub properties: PropertyStore,

    /// Internal index which will be managed by parent molecule
    index: AtomIndex,

    // TODO: remove?
    label: Option<String>,

    /// private atom data independent of the molecule
    data: AtomData,
}

impl Default for Atom {
    fn default() -> Self {
        Atom {
            index: AtomIndex::new(0),
            data: AtomData::default(),
            label: None,

            properties: PropertyStore::new(),
        }
    }
}

impl Atom {
    /// Return element symbol
    pub fn symbol(&self) -> &str {
        self.data.kind.symbol()
    }

    /// Return atomic number
    pub fn number(&self) -> usize {
        self.data.kind.number()
    }

    // FIXME: return &str?
    /// Return element name
    pub fn name(&self) -> String {
        self.data.kind.name()
    }

    /// Provide read-only access to atom index
    pub fn index(&self) -> AtomIndex {
        self.index
    }

    /// Return atom position in 3D Cartesian coordinates
    pub fn position(&self) -> Point3D {
        self.data.position
    }

    /// Set atom position in 3D Cartesian coordinates
    pub fn set_position(&mut self, p: Point3D) {
        self.data.position = p;
    }

    /// Vector quantity equal to the product of mass and velocity.
    pub fn momentum(&mut self) -> Point3D {
        self.data.momentum
    }

    /// TODO: momentum, momenta
    pub fn set_momentum(&mut self, m: Point3D) {
        self.data.momentum = m;
    }

    /// Set atom label
    pub fn set_label(&mut self, lbl: &str) {
        self.label = Some(lbl.into());
    }

    /// Return the user defined atom label, if not return the default (symbol + index, e.g Fe120)
    pub fn label(&self) -> String {
        if let Some(ref l) = self.label {
            return l.to_owned();
        }

        // default atom label: symbol + index
        // example: Fe120
        // counting from 1 instead of 0
        format!("{}{}", self.symbol(), self.index.index() + 1)
    }

    /// Return a list of atoms bonded to current atom
    ///
    /// # Parameters
    ///
    /// parent: the molecule hosting the atom
    ///
    pub fn neighbors<'a>(&self, parent: &'a Molecule) -> Vec<&'a Atom>{
        let atoms: Vec<_> = parent.graph
            .neighbors(self.index)
            .map(|n| &parent.graph[n])
            .collect();

        atoms
    }
}
// 150189fd-57d9-4e19-a888-d64497f5ba7e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b6d1e417-27da-4384-879a-db28960ed161][b6d1e417-27da-4384-879a-db28960ed161]]
use geometry::euclidean_distance;

impl Atom {
    /// Return the distance to other atom
    pub fn distance(&self, other: &Atom) -> f64 {
        euclidean_distance(self.position(), other.position())
    }
}
// b6d1e417-27da-4384-879a-db28960ed161 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::d333cb1f-e622-462f-a892-4906c85b7da0][d333cb1f-e622-462f-a892-4906c85b7da0]]
/// Atom specific data independent of the molecule
#[derive(Debug, Clone)]
pub struct AtomData {
    /// Atom type, could be an element or a pseudo-atom
    kind: AtomKind,
    /// Cartesian coordinates
    position: Point3D,
    /// Atom nick name
    name: String,
    /// Vector quantity equal to the derivative of the position vector with respect to time
    velocity: Point3D,
    /// Atomic mass
    mass: f64,
    /// Atomic momentum vector
    momentum: Point3D,
    /// Atomic partial charge
    partial_charge: f64,

}

impl Default for AtomData {
    fn default()  -> Self {
        AtomData {
            kind: Element(6),   // carbon atom
            position: [0.0; 3],
            momentum: [0.0; 3],
            velocity: [0.0; 3],
            partial_charge: 0.0,
            name: "carbon".into(),
            // FIXME
            mass: 6.0,
        }
    }
}

impl AtomData {
    pub fn new() -> Self {
        AtomData::default()
    }

    /// Set atom position
    pub fn position(&mut self, x: f64, y: f64, z: f64) -> &mut Self {
        self.position = [x, y, z]; self
    }

    /// Set atom kind using element number
    pub fn element(&mut self, n: usize) -> &mut Self {
        self.kind = Element(n); self
    }

    /// Set atom kind using element symbol
    pub fn symbol<T: Into<String>>(&mut self, s: T) -> &mut Self {
        self.kind = atom_kind_from_string(s.into()); self
    }

    pub fn momentum(&mut self, x: f64, y: f64, z: f64) -> &mut Self {
        self.momentum = [x, y, z]; self
    }

    /// return a new `Atom` struct
    pub fn finish(&self) -> Atom {
        let mut atom = Atom::default();
        atom.data = self.clone();

        atom
    }
}

impl Atom {
    pub fn new<T: Into<String>>(s: T, p: Point3D) -> Self {
        AtomData::new()
            .symbol(s)
            .position(p[0], p[1], p[2])
            .finish()
    }

    // FIXME
    pub fn set_symbol<T: Into<String>>(&mut self, symbol: T) {
        self.data.kind = atom_kind_from_string(symbol.into());
    }

    pub fn build() -> AtomData {
        AtomData::new()
    }
}

#[test]
fn test_atom_builder() {
    let a = Atom::build()
        .position(0.0, 0.0, 1.2)
        .symbol("Fe")
        .element(13)
        .momentum(0.2, 0.2, 0.3)
        .finish();

    assert_eq!(13, a.number());
}
// d333cb1f-e622-462f-a892-4906c85b7da0 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::7e463bf4-a6ab-4648-a8a2-4b2023d1c588][7e463bf4-a6ab-4648-a8a2-4b2023d1c588]]
use std::str::FromStr;

impl FromStr for Atom {
    type Err = Error;

    fn from_str(line: &str) -> Result<Self> {
        let parts: Vec<_> = line.split_whitespace().collect();
        if parts.len() != 4 {
            bail!("Incorrect number of data fields: {:?}", line);
        }

        let sym = parts[0];
        let px: f64 = parts[1].parse()?;
        let py: f64 = parts[2].parse()?;
        let pz: f64 = parts[3].parse()?;

        let mut atom = Atom::new(sym, [px, py, pz]);
        atom.set_label(sym);

        Ok(atom)
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:6} {:-12.6} {:-12.6} {:-12.6}",
               self.symbol(),
               self.data.position[0],
               self.data.position[1],
               self.data.position[2]
        )
    }
}
// 7e463bf4-a6ab-4648-a8a2-4b2023d1c588 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b88435fd-d51c-48b8-880c-425b94b905e9][b88435fd-d51c-48b8-880c-425b94b905e9]]
#[test]
fn test_atom_init() {
    let atom = Atom::default();
    let atom = Atom::new("Fe", [9.3; 3]);
    assert_eq!(9.3, atom.position()[0]);
    assert_eq!("Fe", atom.symbol());
    assert_eq!(26, atom.number());

    let atom = Atom::new("dummy", [9.3; 3]);
    assert_eq!("dummy", atom.symbol());
    assert_eq!(0, atom.number());
}
// b88435fd-d51c-48b8-880c-425b94b905e9 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::cfdf0fc1-97a2-4da4-b0bb-a9baee31d275][cfdf0fc1-97a2-4da4-b0bb-a9baee31d275]]
#[test]
fn test_atom_string_conversion() {
    let line = "H 1.0 1.0 1.0";
    let a: Atom = line.parse().unwrap();
    assert_eq!(1, a.number());
    let line = a.to_string();
    let b: Atom = line.parse().unwrap();
    assert_eq!(a.symbol(), b.symbol());
    assert_eq!(a.position(), b.position());
    let line = "24 0.124 1.230 2.349";
    let a: Atom = line.parse().unwrap();
}
// cfdf0fc1-97a2-4da4-b0bb-a9baee31d275 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::7ff70329-69ef-4221-a539-fb097258d0a6][7ff70329-69ef-4221-a539-fb097258d0a6]]
/// https://en.wikipedia.org/wiki/Bond_order
#[derive(Debug, Clone, Copy, PartialOrd, PartialEq)]
pub enum BondKind {
    Dummy,
    Partial,
    Single,
    Aromatic,
    Double,
    Triple,
    Quadruple,
}

/// There is a chemical bond between two atoms or groups of atoms in the case
/// that the forces acting between them are such as to lead to the formation of
/// an aggregate with sufficient stability to make it convenient for the chemist
/// to consider it as an independent 'molecular species'.
///
/// # Reference
/// https://goldbook.iupac.org/html/B/B00697.html
///
#[derive(Debug, Clone)]
pub struct Bond {
    pub kind : BondKind,
    pub name : String,

    /// will be managed by molecule
    index: BondIndex,

    /// set this attribute for arbitrary bond order
    order    : Option<f64>,
}

impl Default for Bond {
    fn default() -> Self {
        Bond {
            order : None,
            kind  : BondKind::Single,

            // private
            name  : String::default(),
            index : BondIndex::new(0),
        }
    }
}

impl Bond {
    pub fn new(order: f64) -> Self {
        debug_assert!(order >= 0.0);

        Bond {
            order: Some(order),
            ..Default::default()
        }
    }

    /// Return bond order
    pub fn order(&self) -> f64 {
        if let Some(order) = self.order {
            order
        } else {
            match self.kind {
                BondKind::Dummy     => 0.0,
                BondKind::Partial   => 0.5,
                BondKind::Single    => 1.0,
                BondKind::Aromatic  => 1.5,
                BondKind::Double    => 2.0,
                BondKind::Triple    => 3.0,
                BondKind::Quadruple => 4.0,
            }
        }
    }

    /// Create a single bond
    pub fn single() -> Self {
        Bond {
            kind: BondKind::Single,
            ..Default::default()
        }
    }

    /// Create a double bond
    pub fn double() -> Self {
        Bond {
            kind: BondKind::Double,
            ..Default::default()
        }
    }

    /// Create a triple bond
    pub fn triple() -> Self {
        Bond {
            kind: BondKind::Triple,
            ..Default::default()
        }
    }

    /// Create an aromatic bond
    pub fn aromatic() -> Self {
        Bond {
            kind: BondKind::Aromatic,
            ..Default::default()
        }
    }

    /// Create a weak bond
    pub fn partial() -> Self {
        Bond {
            kind: BondKind::Partial,
            ..Default::default()
        }
    }

    /// Create a quadruple bond
    pub fn quadruple() -> Self {
        Bond {
            kind: BondKind::Quadruple,
            ..Default::default()
        }
    }

    /// Create a dummy bond
    pub fn dummy() -> Self {
        Bond {
            kind: BondKind::Dummy,
            ..Default::default()
        }
    }

    /// Read-only access of bond index
    pub fn index(&self) -> BondIndex {
        self.index
    }

    /// Return indices to atoms in parent molecule hold by the bond
    pub fn partners<'a>(&self, parent: &'a Molecule) -> Option<(&'a Atom, &'a Atom)> {
        if let Some((n1, n2)) = parent.graph.edge_endpoints(self.index) {
            let a1 = &parent.graph[n1];
            let a2 = &parent.graph[n2];
            return Some((a1, a2));
        }

        None
    }
}
// 7ff70329-69ef-4221-a539-fb097258d0a6 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::486bd5a4-e762-46bf-a237-e692393a795d][486bd5a4-e762-46bf-a237-e692393a795d]]
#[test]
fn test_bond() {
    let b = Bond::default();
    let b = Bond::new(1.5);
    assert_eq!(1.5, b.order());
}
// 486bd5a4-e762-46bf-a237-e692393a795d ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::942dedaa-9351-426e-9be9-cdb640ec2b75][942dedaa-9351-426e-9be9-cdb640ec2b75]]
pub type MolGraph = StableUnGraph<Atom, Bond>;
pub type AtomIndex = NodeIndex;
pub type BondIndex = EdgeIndex;

/// Molecule is the most important data structure in gchemol, which repsents
/// "any singular entity, irrespective of its nature, used to concisely express
/// any type of chemical particle that can exemplify some process: for example,
/// atoms, molecules, ions, etc. can all undergo a chemical reaction". Molecule
/// may have chemical bonds between atoms.
///
/// Reference
/// ---------
/// 1. http://goldbook.iupac.org/M03986.html
/// 2. https://en.wikipedia.org/wiki/Molecular_entity
///
#[derive(Debug, Clone)]
pub struct Molecule {
    /// Molecule name
    pub name: String,
    /// core data in graph
    pub graph: MolGraph,
    /// Crystalline lattice for structure using periodic boundary conditions
    pub lattice: Option<Lattice>,
    /// Arbitrary molecular property stored in key-value pair. Key is a string
    /// type, but it is the responsibility of the setter/getter to interpret the
    /// value.
    pub properties: PropertyStore,
}

impl Default for Molecule {
    fn default() -> Self {
        let graph = MolGraph::default();
        Molecule {
            name: "default".to_string(),
            graph: graph,
            lattice: None,
            properties: PropertyStore::new(),
        }
    }
}

impl Molecule {
    /// Create a new empty molecule with specific name
    pub fn new(name: &str) -> Self {
        Molecule {
            name: name.to_string(),
            ..Default::default()
        }
    }

    /// Return the number of atoms in the molecule.
    pub fn natoms(&self) -> usize {
        self.graph.node_count()
    }

    /// Return the number of bonds in the molecule.
    pub fn nbonds(&self) -> usize {
        self.graph.edge_count()
    }

    /// Construct from an existing graph
    pub fn from_graph(graph: MolGraph) -> Self{
        Molecule {
            graph: graph,
            ..Default::default()
        }
    }

    /// A convenient alias of molecular name
    pub fn title(&self) -> String {
        self.name.to_owned()
    }

    /// Return the sites (AtomIndex) that hosting atoms.
    pub fn sites(&self) -> Vec<AtomIndex> {
        self.graph.node_indices().collect()
    }

    /// Return an iterator over the atoms in the molecule.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.graph.node_indices().map(move |n| &self.graph[n])
    }

    /// Return an iterator over the bonds in the molecule.
    pub fn bonds(&self) -> impl Iterator<Item = &Bond> {
        self.graph.edge_indices().map(move |e| &self.graph[e])
    }

    /// Return positions of all atoms in the molecule.
    pub fn positions(&self) -> Vec<Point3D> {
        self.atoms().map(|ref a| a.position()).collect()
    }

    /// Return symbols of all atoms in the molecule.
    pub fn symbols(&self) -> Vec<&str> {
        self.atoms().map(|ref a| a.symbol()).collect()
    }

    /// Return element numbers of all atoms in the molecule.
    pub fn numbers(&self) -> Vec<usize> {
        self.atoms().map(|ref a| a.number()).collect()
    }

    // FIXME: opt performance
    /// Set positions of atoms
    pub fn set_positions(&mut self, positions: Points) -> Result<()>
    {
        let indices = self.sites();
        if indices.len() != positions.len() {
            bail!("the number of cartesian coordinates is different from the number of atoms in molecule.")
        }

        for (&index, position) in indices.iter().zip(positions) {
            let mut atom = &mut self.graph[index];
            atom.set_position(position);
        }

        Ok(())
    }

    // FIXME: add a test
    /// Set element symbols
    pub fn set_symbols<'a, I>(&mut self, symbols: I) -> Result<()>
    where
        I: IntoIterator,
        I::Item: Into<String>,
    {
        let indices = self.sites();
        let symbols = symbols.into_iter();

        for (&index, symbol) in indices.iter().zip(symbols) {
            let mut atom = &mut self.graph[index];
            atom.set_symbol(symbol);
        }

        Ok(())
    }
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::be29e151-18c6-43cb-9586-aba0e708d38c][be29e151-18c6-43cb-9586-aba0e708d38c]]
pub trait IntoAtomIndex {
    fn into_atom_index(&self) -> AtomIndex;
}

impl IntoAtomIndex for usize {
    fn into_atom_index(&self) -> AtomIndex {
        // atom index counting from zero
        AtomIndex::new(*self)
    }
}

impl IntoAtomIndex for AtomIndex {
    fn into_atom_index(&self) -> AtomIndex {
        *self
    }
}

pub trait IntoBondIndex {
    fn into_bond_index(&self) -> BondIndex;
}

impl IntoBondIndex for usize {
    fn into_bond_index(&self) -> BondIndex {
        // atom index counting from zero
        BondIndex::new(*self)
    }
}

impl IntoBondIndex for BondIndex {
    fn into_bond_index(&self) -> BondIndex {
        *self
    }
}
// be29e151-18c6-43cb-9586-aba0e708d38c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::80dcc47b-b7dc-4ba7-a9d6-a567831bae93][80dcc47b-b7dc-4ba7-a9d6-a567831bae93]]
impl Molecule {
    // TODO
    /// Return molecule net charge
    pub fn charge(&self) -> usize {
        unimplemented!()
    }

    // TODO
    /// Return spin multiplicity of the molecule
    pub fn spin_multiplicity(&self) -> usize {
        unimplemented!()
    }

    // TODO
    /// Return the number of electrons in the system, based on the atomic numbers and
    /// molecular formal charge
    pub fn nelectrons(&self) -> usize {
        unimplemented!()
    }
}
// 80dcc47b-b7dc-4ba7-a9d6-a567831bae93 ends here
