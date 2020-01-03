// header

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*header][header:1]]
//===============================================================================#
//   DESCRIPTION:  molecule object repsented in graph data structure
//
//       OPTIONS:  ---
//  REQUIREMENTS:  ---
//         NOTES:  based on petgraph
//        AUTHOR:  Wenping Guo <ybyygu@gmail.com>
//       LICENCE:  GPL version 3
//       CREATED:  <2018-04-12 Thu 15:48>
//       UPDATED:  <2020-01-03 Fri 12:57>
//===============================================================================#
// header:1 ends here

// base

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*base][base:1]]
use std::collections::HashMap;
use std::convert;

use petgraph;
use petgraph::prelude::*;
use serde::{Deserialize, Serialize};

use self::property::PropertyStore;
use crate::core_utils::*;
use crate::lattice::Lattice;

mod clean;
mod connect;
mod edit;
mod formula;
mod fragment;
mod order;
mod property;
mod view;

#[cfg(test)]
mod test;

pub(crate) type Point3D = [f64; 3];
pub(crate) type Points = Vec<Point3D>;
// base:1 ends here

// globals

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*globals][globals:1]]
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
    ("Uuo", "ununoctium"),
];
// globals:1 ends here

// base

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*base][base:1]]
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub enum AtomKind {
    /// physical elements
    Element(usize),
    /// a dummy-atom is not a real atom
    Dummy(String),
}

impl AtomKind {
    pub fn symbol(&self) -> &str {
        match self {
            &Element(num) => ELEMENT_DATA[num - 1].0,
            &Dummy(ref sym) => sym,
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
            &Element(num) => String::from(ELEMENT_DATA[num - 1].1),
            &Dummy(ref sym) => format!("dummy atom {}", sym),
        }
    }
}

impl fmt::Display for AtomKind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            &Element(num) => write!(f, "{:}", self.symbol()),
            &Dummy(ref sym) => write!(f, "{:}", self.symbol()),
        }
    }
}
// base:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*test][test:1]]
use self::AtomKind::{Dummy, Element};

/// Return AtomKind using common sense
pub fn atom_kind_from_string<T: Into<String>>(sym: T) -> AtomKind {
    let sym = sym.into();

    // element specified in number
    if let Ok(x) = sym.parse::<usize>() {
        return Element(x);
    }

    // element specified in symbol or long name
    let _sym = sym.to_uppercase();
    for (i, &(s, n)) in ELEMENT_DATA.iter().enumerate() {
        if s.to_uppercase() == _sym || n.to_uppercase() == _sym {
            return Element(i + 1);
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

    let k = atom_kind_from_string("SI");
    assert_eq!(k.symbol(), "Si");
}
// test:1 ends here

// base

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*base][base:1]]
use std::hash::Hash;
use std::cmp::Ordering;

/// Atom is the smallest particle still characterizing a chemical element.
///
/// # Reference
///
/// https://goldbook.iupac.org/html/A/A00493.html
///
#[derive (Debug, Clone, Deserialize, Serialize)]
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
    pub fn set_position(&mut self, p: [f64; 3]) {
        self.data.position = p;
    }

    /// Vector quantity equal to the product of mass and velocity.
    pub fn momentum(&self) -> Point3D {
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
// base:1 ends here

// geometry

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*geometry][geometry:1]]
use crate::geometry::prelude::euclidean_distance;

impl Atom {
    /// Return the distance to other atom
    pub fn distance(&self, other: &Atom) -> f64 {
        euclidean_distance(self.position(), other.position())
    }
}
// geometry:1 ends here

// atom data builder
// - 目的: 方便的设置不常用的Atom属性, 固定外围参数的使用接口, 隔离Atom内部属性的实
//   现细节.
// - 使用builder pattern

// References
// - servo: [[https://doc.servo.org/src/cookie/builder.rs.html#35-38][builder.rs.html -- source]]


// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*atom data builder][atom data builder:1]]
/// Atom specific data independent of the molecule
#[derive(Debug, Clone, Deserialize, Serialize)]
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
    fn default() -> Self {
        AtomData {
            kind: Element(6), // carbon atom
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
        self.position = [x, y, z];
        self
    }

    /// Set atom kind using element number
    pub fn element(&mut self, n: usize) -> &mut Self {
        self.kind = Element(n);
        self
    }

    /// Set atom kind using element symbol
    pub fn symbol<T: Into<String>>(&mut self, s: T) -> &mut Self {
        self.kind = atom_kind_from_string(s.into());
        self
    }

    pub fn momentum(&mut self, x: f64, y: f64, z: f64) -> &mut Self {
        self.momentum = [x, y, z];
        self
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

    // pub fn set_symbol<T: Into<String>>(&mut self, symbol: T) {
    //     self.data.kind = atom_kind_from_string(symbol.into());
    // }

    pub fn set_symbol(&mut self, symbol: &str) {
        self.data.kind = atom_kind_from_string(symbol);
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
// atom data builder:1 ends here

// convert

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*convert][convert:1]]
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
        write!(
            f,
            "{:6} {:-12.6} {:-12.6} {:-12.6}",
            self.symbol(),
            self.data.position[0],
            self.data.position[1],
            self.data.position[2]
        )
    }
}

impl<S, C> std::convert::From<(S, C)> for Atom
where
    S: Into<String>,
    C: Into<[f64; 3]>,
{
    fn from(p: (S, C)) -> Self {
        let s = p.0.into();
        let c = p.1.into();
        Atom::new(s, c)
    }
}
// convert:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*test][test:1]]
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
// test:1 ends here

// src

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*src][src:1]]
/// https://en.wikipedia.org/wiki/Bond_order
#[derive(Debug, Clone, Copy, PartialOrd, PartialEq, Deserialize, Serialize)]
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
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Bond {
    pub kind: BondKind,
    pub name: String,

    /// will be managed by molecule
    index: BondIndex,

    /// set this attribute for arbitrary bond order
    order: Option<f64>,
}

impl Default for Bond {
    fn default() -> Self {
        Bond {
            order: None,
            kind: BondKind::Single,

            // private
            name: String::default(),
            index: BondIndex::new(0),
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
                BondKind::Dummy => 0.0,
                BondKind::Partial => 0.5,
                BondKind::Single => 1.0,
                BondKind::Aromatic => 1.5,
                BondKind::Double => 2.0,
                BondKind::Triple => 3.0,
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
// src:1 ends here

// test

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*test][test:1]]
#[test]
fn test_bond() {
    let b = Bond::default();
    let b = Bond::new(1.5);
    assert_eq!(1.5, b.order());
}
// test:1 ends here

// molecule

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*molecule][molecule:1]]
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
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Molecule {
    /// Molecular name. The default value is `untitled`. This field may be used
    /// to store a registry number or other identifier, instead of a common
    /// name.
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
            name: String::new(),
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
    pub fn from_graph(graph: MolGraph) -> Self {
        Molecule {
            graph: graph,
            ..Default::default()
        }
    }

    /// Return the name of the molecule, while is typpically modified for safely
    /// storing in various molecular file formats.
    pub fn title(&self) -> String {
        let tlines: Vec<_> = self.name.lines().collect();
        if tlines.is_empty() {
            "untitled".to_owned()
        } else {
            tlines[0].trim().to_owned()
        }
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
}
// molecule:1 ends here

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*molecule][molecule:2]]
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
// molecule:2 ends here

// TODO charge, spin, magmon
// 定义分子体系的电荷, 自旋等参数. ase里每个原子都有一个magmon的参数.


// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*charge, spin, magmon][charge, spin, magmon:1]]
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
// charge, spin, magmon:1 ends here
