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
//       UPDATED:  <2018-06-21 Thu 16:38>
//===============================================================================#

use std::collections::HashMap;
use petgraph;
use petgraph::prelude::*;
use std::convert;
use quicli::prelude::*;

use {
    Point3D,
    Points,
    lattice::Lattice,
};
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

/// Molecule is the most important data structure in gchemol, which repsents one or
/// more atoms held together by chemical bonds.
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

    /// Mapping user defined atom index to internal graph node index
    atom_indices: HashMap<String, NodeIndex>,
    /// Mapping bond tuple to EdgeIndex
    bond_indices: HashMap<[String; 2], EdgeIndex>,
    /// User defined atom labels
    pub atom_labels: HashMap<AtomIndex, String>,
}

impl Default for Molecule {
    fn default() -> Self {
        let graph = MolGraph::default();
        Molecule {
            name: "default".to_string(),
            graph: graph,
            lattice: None,
            properties: PropertyStore::new(),
            atom_indices: HashMap::new(),
            bond_indices: HashMap::new(),
            atom_labels: HashMap::new(),
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

    /// Return symbols of all  atoms in the molecule.
    pub fn symbols(&self) -> Vec<&str> {
        self.atoms().map(|ref a| a.symbol()).collect()
    }

    // FIXME: opt performance
    /// Set positions of atoms
    pub fn set_positions(&mut self, positions: Points) -> Result<()>
    {
        let indices: Vec<_> = self.graph.node_indices().collect();
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
        let indices: Vec<_> = self.graph.node_indices().collect();
        let symbols = symbols.into_iter();

        for (&index, symbol) in indices.iter().zip(symbols) {
            let mut atom = &mut self.graph[index];
            atom.set_symbol(symbol);
        }

        Ok(())
    }
}
// 942dedaa-9351-426e-9be9-cdb640ec2b75 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ab071843-5a26-4f16-9068-17da002d5a10][ab071843-5a26-4f16-9068-17da002d5a10]]
impl Molecule {
    /// Return fractional coordinates relative to unit cell.
    pub fn scaled_positions(&self) -> Option<Vec<Point3D>> {
        if let Some(mut lat) = self.lattice {
            let mut fxyzs = vec![];
            for a in self.atoms() {
                let xyz = a.position();
                let fxyz = lat.to_frac(xyz);
                fxyzs.push(fxyz)
            }
            Some(fxyzs)
        } else {
            None
        }
    }

    /// Set positions relative to unit cell.
    pub fn set_scaled_positions(&mut self, scaled: &Points) -> Result<()> {
        if let Some(mut lat) = self.lattice {
            let mut positions = vec![];
            for &p in scaled {
                let xyz = lat.to_cart(p);
                positions.push(p);
            }

            self.set_positions(positions)
        } else {
            bail!("cannot set scaled positions for aperiodic structure")
        }
    }
}
// ab071843-5a26-4f16-9068-17da002d5a10 ends here

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

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a2f1dbf4-c41d-4ca4-986e-5e8cfb3d5d08][a2f1dbf4-c41d-4ca4-986e-5e8cfb3d5d08]]
use serde::{
    de::DeserializeOwned,
    ser::Serialize,
};

use serde_json;
use std::result;

/// A container storing extra information managed as key/value pairs
#[derive(Debug, Clone)]
pub struct PropertyStore {
    data: HashMap<String, String>,
}

impl PropertyStore {
    fn new() -> Self {
        PropertyStore {
            data: HashMap::new(),
        }
    }

    /// retrieve property associated with the `key`
    pub fn load<D: DeserializeOwned>(&self, key: &str) -> result::Result<D, serde_json::Error> {
        let serialized = self.data.get(key).unwrap();
        serde_json::from_str(&serialized)
    }

    /// store property associatd with a `key`
    pub fn store<D: Serialize>(&mut self, key: &str, value: D) {
        let serialized = serde_json::to_string(&value).unwrap();
        self.data.insert(key.into(), serialized);
    }

    pub fn discard(&mut self, key: &str) {
        self.data.remove(key.into());
    }
}

#[test]
fn test_atom_store() {
    let mut x = PropertyStore::new();
    let d = [1, 2, 3];
    x.store("k", d);
    let x: [usize; 3] = x.load("k").unwrap();
    assert_eq!(d, x);
}
// a2f1dbf4-c41d-4ca4-986e-5e8cfb3d5d08 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::9924e323-dd02-49d0-ab07-41208114546f][9924e323-dd02-49d0-ab07-41208114546f]]
impl Molecule {
    /// Add a single atom into molecule
    /// Return an index to atom counting from 1
    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        let n = self.graph.add_node(atom);
        let atom = self.get_atom_mut(n).unwrap();
        atom.index = n;

        n
    }

    /// Remove an atom from the molecule.
    /// Return the removed atom if it exists, or return None.
    pub fn remove_atom<T: IntoAtomIndex>(&mut self, index: T) -> Option<Atom> {
        let n = index.into_atom_index();
        self.graph.remove_node(n)
    }

    /// access atom by atom index
    pub fn get_atom<T: IntoAtomIndex>(&self, index: T) -> Option<&Atom> {
        let n = index.into_atom_index();
        self.graph.node_weight(n)
    }

    /// mutable access to atom by atom index
    pub fn get_atom_mut<T: IntoAtomIndex>(&mut self, index: T) -> Option<&mut Atom> {
        let n = index.into_atom_index();
        self.graph.node_weight_mut(n)
    }

    /// Add a single bond into molecule specified by atom indices
    /// Will panic if corresponding atoms does not exist
    /// The existing bond data will be replaced if n1 already bonded with n2
    /// Return a bond index pointing to bond data
    pub fn add_bond(&mut self, n1: AtomIndex, n2: AtomIndex, bond: Bond) -> BondIndex {
        let e = self.graph.update_edge(n1, n2, bond);

        // cache the pair of atoms
        let mut bond = &mut self.graph[e];
        bond.index = e;

        e
    }

    /// Access bond by bond index
    pub fn get_bond<T: IntoBondIndex>(&self, e: T) -> Option<&Bond> {
        let e = e.into_bond_index();
        self.graph.edge_weight(e)
    }

    /// Get the bond index between two atoms.
    /// Return None if not found
    fn bond_index_between(&self, n1: AtomIndex, n2: AtomIndex) -> Option<BondIndex> {
        self.graph.find_edge(n1, n2)
    }

    /// Return any bond bween two atoms.
    /// Return None if it does not exist
    pub fn get_bond_between<T: IntoAtomIndex>(&self, index1: T, index2: T) -> Option<&Bond> {
        let n1 = index1.into_atom_index();
        let n2 = index2.into_atom_index();

        if let Some(e) = self.bond_index_between(n1, n2) {
            self.get_bond(e)
        } else {
            None
        }
    }

    /// Remove a bond specified by bond index
    /// Return the removed bond if it exists, or return None
    pub fn remove_bond<T: IntoBondIndex>(&mut self, b: T) -> Option<Bond>{
        let b = b.into_bond_index();
        self.graph.remove_edge(b)
    }

    /// Remove any bond bween two atoms
    /// Return the removed bond if it exists, or return None
    pub fn remove_bond_between<T: IntoAtomIndex>(&mut self, index1: T, index2: T) -> Option<Bond> {
        let n1 = index1.into_atom_index();
        let n2 = index2.into_atom_index();

        if let Some(e) = self.bond_index_between(n1, n2) {
            self.remove_bond(e)
        } else {
            None
        }
    }

    /// Removes all bonds between two selections to respect pymol's unbond command.
    ///
    /// Parameters
    /// ----------
    /// atom_indices1: the first collection of atoms
    ///
    /// atom_indices2: the other collection of atoms
    ///
    /// Reference
    /// ---------
    /// https://pymolwiki.org/index.php/Unbond
    ///
    pub fn unbond(&mut self, atom_indices1: Vec<AtomIndex>, atom_indices2: Vec<AtomIndex>) {
        for &index1 in atom_indices1.iter() {
            for &index2 in atom_indices2.iter() {
                self.remove_bond_between(index1, index2);
            }
        }
    }

    /// mutable access to bond by bond index
    pub fn get_bond_mut<T: IntoBondIndex>(&mut self, b: T) -> Option<&mut Bond> {
        let b = b.into_bond_index();
        self.graph.edge_weight_mut(b)
    }

    /// Update atom indices
    pub fn reorder(&mut self) {
        let ns: Vec<_> = self.graph.node_indices().collect();
        for (i, &n) in ns.iter().enumerate() {
            let atom = &mut self.graph[n];
            atom.index = n;
        }
    }
}
// 9924e323-dd02-49d0-ab07-41208114546f ends here

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

// [[file:~/Workspace/Programming/gchemol/gchemol.note::72dd0c31-26e5-430b-9f67-1c5bd5220a84][72dd0c31-26e5-430b-9f67-1c5bd5220a84]]
impl Molecule {
    /// add many atoms from a hashmap
    pub fn add_atoms_from(&mut self, atoms: HashMap<&str, Atom>) -> Result<()>{
        for (k, a) in atoms {
            let n = self.add_atom(a);
            self.atom_labels.insert(n, k.into());
            self.atom_indices.insert(k.into(), n);
        }

        Ok(())
    }

    pub fn add_bonds_from(&mut self, bonds: HashMap<(String, String), Bond>) -> Result<()>{
        for ((ki, kj), b) in bonds {
            if ki == kj {
                bail!("Bonding with self is not allowed.");
            }

            if let Some(&ni) = self.atom_indices.get(&ki) {
                if let Some(&nj) = self.atom_indices.get(&kj) {
                    // get bond key from the user
                    let kij = if ki < kj {[ki, kj]} else {[kj, ki]};
                    // get internal graph edge index
                    let e = self.add_bond(ni, nj, b);
                    self.bond_indices.insert(kij, e);
                } else {
                    bail!("atom {} is not presented in molecule.", kj);
                }
            } else {
                bail!("atom {} is not presented in molecule.", ki);
            }
        }

        Ok(())
    }

    pub fn remove_atoms_from(&mut self) {
        unimplemented!()
    }

    pub fn remove_bonds_from(&mut self) {
        unimplemented!()
    }
}
// 72dd0c31-26e5-430b-9f67-1c5bd5220a84 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::5754ca07-a93d-47e5-8256-c7236777b2ee][5754ca07-a93d-47e5-8256-c7236777b2ee]]
impl Molecule {
    /// Set periodic lattice
    pub fn set_lattice(&mut self, lat: Lattice) {
        self.lattice = Some(lat);
    }

    /// Unbuild current crystal structure leaving a nonperiodic structure
    pub fn unbuild_crystal(&mut self) {
        self.lattice = None
    }
}
// 5754ca07-a93d-47e5-8256-c7236777b2ee ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::1f8f3a9a-dc8d-4a1e-816b-c5228ddf0607][1f8f3a9a-dc8d-4a1e-816b-c5228ddf0607]]
impl Molecule {
    // FIXME: if PBC
    pub fn distance_matrix(&self) -> Vec<Vec<f64>>{
        let positions = self.positions();
        get_distance_matrix(positions)
    }

    // TODO: improve performance
    /// Return the distance between `atom i` and `atom j`.
    ///
    /// Force periodic structure, this method will return the distance under the
    /// minimum image convention.
    pub fn distance(&self, i: AtomIndex, j: AtomIndex) -> Option<f64> {
        if let Some(ai) = self.get_atom(i) {
            if let Some(aj) = self.get_atom(j) {
                if let Some(mut lat) = self.lattice {
                    let pi = ai.position();
                    let pj = aj.position();
                    let dij = lat.distance(pi, pj);
                    Some(dij)
                } else {
                    Some(ai.distance(aj))
                }
            } else {
                None
            }
        } else {
            None
        }
    }
}
// 1f8f3a9a-dc8d-4a1e-816b-c5228ddf0607 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3][a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3]]
use geometry::get_distance_matrix;
use data::guess_bond_kind;

impl Molecule {
    /// Access the bonded atom indices for bond `b`
    pub fn partners<T: IntoBondIndex>(&self, b: &T) -> Option<(AtomIndex, AtomIndex)>{
        let b = b.into_bond_index();
        self.graph.edge_endpoints(b)
    }

    /// Return all connected atoms with `a`
    pub fn neighbors(&self, a: AtomIndex) -> Vec<AtomIndex> {
        self.graph.neighbors(a).collect()
    }

    /// Removes all existing bonds between atoms
    pub fn unbound(&mut self) {
        self.graph.clear_edges();
    }

    /// Redefine bonds from distances based on predefined bonding lengths
    pub fn rebond(&mut self) {
        let n = self.natoms();
        let mut pairs = vec![];
        for ni in self.graph.node_indices() {
            for nj in self.graph.node_indices() {
                if ni < nj {
                    let ai = &self.graph[ni];
                    let aj = &self.graph[nj];
                    let bk = guess_bond_kind(ai, aj);
                    if bk != BondKind::Dummy {
                        let mut bond = Bond::default();
                        bond.kind = bk;
                        pairs.push((ni, nj, bond));
                    }
                }
            }
        }

        for (i, j, b) in pairs {
            self.add_bond(i, j, b);
        }
    }
}
// a1ee57e8-ac54-4e78-9e8a-a5b5bf11f0e3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ec7b11d2-6f13-49fd-b253-af4b213b49a3][ec7b11d2-6f13-49fd-b253-af4b213b49a3]]
impl Molecule {
    /// Return an iterator of all atoms connected to a.
    fn connected() {
        //
    }
}
// ec7b11d2-6f13-49fd-b253-af4b213b49a3 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::2a27ca30-0a99-4d5d-b544-5f5900304bbb][2a27ca30-0a99-4d5d-b544-5f5900304bbb]]
use petgraph::algo;
use rand::{thread_rng, Rng};

const EPSILON: f64 = 1.0E-6;

impl Molecule {
    /// Set atom position
    pub fn set_position(&mut self, index: AtomIndex, position: Point3D) {
        let atom = &mut self.graph[index];
        atom.set_position(position);
    }

    /// Return the shortest distance numbered in bonds between two atoms
    /// Return None if them are not connected
    pub fn nbonds_between(&self, index1: AtomIndex, index2: AtomIndex) -> Option<usize> {
        let path = algo::astar(&self.graph,
                               index1,
                               |finish| finish == index2,
                               |e| 1,
                               |_| 0);
        if let Some((n, _)) = path {
            Some(n)
        } else {
            None
        }
    }

    /// Translate atomic positions by a displacement
    pub fn translate(&mut self, displacement: Point3D) {
        let nodes: Vec<_> = self.graph.node_indices().collect();
        for n in nodes {
            let mut atom = &mut self.graph[n];
            let mut position = atom.position();
            for v in 0..3 {
                position[v] += displacement[v];
            }
            atom.set_position(position);
        }
    }

    /// Return the center of mass of molecule (COM).
    pub fn center_of_mass(&self) -> Point3D {
        unimplemented!()
    }

    /// Return the center of geometry of molecule (COG).
    pub fn center_of_geometry(&self) -> Point3D {
        let mut p = [0.0; 3];
        for [x, y, z] in self.positions() {
            p[0] += x;
            p[1] += y;
            p[2] += z;
        }

        let n = self.natoms() as f64;
        p[0] /= n;
        p[1] /= n;
        p[2] /= n;

        p
    }

    /// Center the molecule around its center of geometry
    pub fn recenter(&mut self) {
        let mut p = self.center_of_geometry();
        for i in 0..3 {
            p[i] *= -1.0;
        }
        self.translate(p);
    }
}

#[test]
fn test_molecule_center() {
    let mol = Molecule::from_file("tests/files/mol2/alanine-gv.mol2").expect("mol2 gv");
    let pc = mol.center_of_geometry();
    let pe = [-2.31413333, -1.24455833,  0.41005833];

    for i in 0..3 {
        assert_relative_eq!(pe[i], pc[i], epsilon=1e-4);
    }

    // TODO
    // center of mass
    // let pe = [-1.8295483 , -1.10700382,  0.25332597];
}
// 2a27ca30-0a99-4d5d-b544-5f5900304bbb ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::e1d0c51a-0dd7-4977-ae54-7928ee46d373][e1d0c51a-0dd7-4977-ae54-7928ee46d373]]
use std::ops::Index;

/// A list-like object providing a convenient view on atoms in molecule
#[derive(Debug, Clone)]
pub struct AtomsView<'a> {
    /// mapping a positive integer to internal graph node index
    mapping: HashMap<usize, AtomIndex>,
    /// parent molecule struct
    parent: &'a Molecule,

    // current position in iteration
    cur: usize,
}

impl<'a> AtomsView<'a> {
    pub fn new(mol: &'a Molecule) -> Self {
        // use a hash map to cache graph node indices
        let mut mapping = HashMap::new();
        let mut i = 1;
        for index in mol.graph.node_indices() {
            mapping.insert(i, index);
            i += 1;
        }

        AtomsView {
            mapping,
            parent: mol,
            cur: 0,
        }
    }
}

/// Index the atoms in `Molecule` by index counting from 1
/// Will panic if index is invalid
impl<'a> Index<usize> for AtomsView<'a>
{
    type Output = Atom;

    fn index(&self, index: usize) -> &Atom {
        let n = self.mapping[&index];
        &self.parent.graph[n]
    }
}

impl<'a> Iterator for AtomsView<'a> {
    type Item = (usize, &'a Atom);

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur >= self.mapping.len() {
            None
        } else {
            self.cur += 1;
            let n = self.mapping[&self.cur];
            let a = &self.parent.graph[n];

            Some((self.cur, &a))
        }
    }
}

#[test]
fn test_atom_view() {
    let mut mol = Molecule::default();
    mol.add_atom(Atom::new("Fe", [0.0; 3]));
    mol.add_atom(Atom::new("C", [0.0; 3]));

    let av = AtomsView::new(&mol);
    assert_eq!("Fe", av[1].symbol());

    // iterate with a index (counting from 1) and an atom object
    for (i, a) in av {
        //
    }
}

impl Molecule {
    pub fn view_atoms(&self) -> AtomsView {
        AtomsView::new(&self)
    }
}
// e1d0c51a-0dd7-4977-ae54-7928ee46d373 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::5916eec2-ec7e-4525-bc6c-fade1d250a16][5916eec2-ec7e-4525-bc6c-fade1d250a16]]
use indexmap::IndexMap;

/// A list-like object providing a convenient view on bonds in molecule
#[derive(Debug, Clone)]
pub struct BondsView<'a> {
    mapping: IndexMap<(usize, usize), BondIndex>,
    parent: &'a Molecule,

    // current position in iteration
    cur: usize,
}

impl<'a> BondsView<'a> {
    pub fn new(mol: &'a Molecule) -> Self {
        // reverse mapping atom id and internal graph node index
        let mut d = HashMap::new();
        let mut i = 1;
        for n in mol.graph.node_indices() {
            d.insert(n, i);
            i += 1;
        }
        // use a hash map to cache graph edge indices
        let mut mapping = indexmap!{};
        for e in mol.graph.edge_indices() {
            let (ni, nj) = mol.graph.edge_endpoints(e).expect("bondview endpoints");
            let ai = d[&ni];
            let aj = d[&nj];
            // make sure ai is always smaller than aj
            if ai < aj {
                mapping.insert((ai, aj), e);
            } else {
                mapping.insert((aj, ai), e);
            }
        }

        BondsView {
            mapping,
            parent: mol,
            cur: 0,
        }
    }
}

impl<'a> Index<(usize, usize)> for BondsView<'a>
{
    type Output = Bond;

    fn index(&self, b: (usize, usize)) -> &Bond {
        // make sure the first index number is always smaller
        let e = if b.0 < b.1 {
            self.mapping[&b]
        } else {
            self.mapping[&(b.1, b.0)]
        };

        &self.parent.graph[e]
    }
}

impl<'a> Iterator for BondsView<'a> {
    type Item = (usize, usize, &'a Bond);

    /// return a tuple in (index_i, index_j, bond)
    /// index_i is always smaller than index_j
    fn next(&mut self) -> Option<Self::Item> {
        if self.cur >= self.mapping.len() {
            None
        } else {
            let (&(i, j), &e) = self.mapping.get_index(self.cur).expect("bondview: get bond by index");
            let b = &self.parent.get_bond(e).expect("bondview: get bond by edge index");
            self.cur += 1;

            if i < j {
                Some((i, j, &b))
            } else {
                Some((j, i, &b))
            }
        }
    }
}

#[test]
fn test_bonds_view() {
    let mut mol = Molecule::default();
    let a1 = mol.add_atom(Atom::new("C", [0.0; 3]));
    let a2 = mol.add_atom(Atom::new("H", [1.0; 3]));
    let a3 = mol.add_atom(Atom::new("H", [2.0; 3]));
    mol.add_bond(a1, a2, Bond::default());
    mol.add_bond(a1, a3, Bond::default());
    // update bond type
    mol.add_bond(a3, a1, Bond::double());
    let bv = BondsView::new(&mol);
    {
        let b12 = &bv[(1, 2)];
        let b13 = &bv[(1, 3)];
    }

    // a list of tuple: (1, 2, Bond)
    let ijbs: Vec<_> = bv.collect();
    assert_eq!(2, ijbs.len());
    let (i, j, _) = ijbs[0];
    assert_eq!((1, 2), (i, j));
}

impl Molecule {
    pub fn view_bonds(&self) -> BondsView {
        BondsView::new(&self)
    }
}
// 5916eec2-ec7e-4525-bc6c-fade1d250a16 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::e2130a32-e39f-4b7b-9014-515f18ff5f48][e2130a32-e39f-4b7b-9014-515f18ff5f48]]
impl Molecule {
    /// Interpolate between this structure and end_structure. Useful for
    /// construction of NEB inputs.
    fn interpolate(&self, end_mol: &Molecule, nimages: usize) -> Vec<Molecule> {
        unimplemented!()
    }
}
// e2130a32-e39f-4b7b-9014-515f18ff5f48 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::82294367-1b69-4638-a70b-fd8daf02ff3e][82294367-1b69-4638-a70b-fd8daf02ff3e]]
type Bounds = HashMap<(AtomIndex, AtomIndex), f64>;

// return distance bounds between atoms
// upper-tri for upper bounds
// lower-tri for lower bounds
fn get_distance_bounds_v1(mol: &Molecule) -> Bounds {
    let mut dm = mol.distance_matrix();
    // max distance between two atoms
    let max_rij = 90.0;

    let mut bounds = HashMap::new();
    let node_indices: Vec<_> = mol.graph.node_indices().collect();
    let nnodes = node_indices.len();
    for node_i in mol.graph.node_indices() {
        for node_j in mol.graph.node_indices() {
            if node_i >= node_j {
                continue;
            }

            let atom_i = &mol.graph[node_i];
            let atom_j = &mol.graph[node_j];

            // use vdw radii as the lower bound for non-bonded pair
            let vri = atom_i.vdw_radius().unwrap();
            let vrj = atom_j.vdw_radius().unwrap();
            let vrij = vri + vrj;

            // use covalent radii as the lower bound for bonded pair
            let cri = atom_i.cov_radius().unwrap();
            let crj = atom_j.cov_radius().unwrap();
            let crij = cri + crj;

            let lij = crij * 0.8;
            let uij = vrij * 0.8;
            let uij = if uij > crij*1.2 {uij} else {crij*1.2};

            // make sure vdw radius larger than covalent radius (usually it is)
            let mut bound = [lij, uij];
            if crij > vrij {
                bound.swap(0, 1);
            }

            let dij = atom_i.distance(atom_j);
            // if i and j is directly bonded
            // set covalent radius as the lower bound
            // or set vdw radius as the lower bound if not bonded
            if let Some(nb) = mol.nbonds_between(node_i, node_j) {
                if nb == 1 {
                    if dij > bound[1] && dij < max_rij {
                        bounds.insert((node_i, node_j), dij);
                        bounds.insert((node_j, node_i), dij);
                    } else {
                        bounds.insert((node_i, node_j), bound[0]);
                        bounds.insert((node_j, node_i), bound[0]);
                    }
                } else if nb == 2 {
                    if dij > bound[1] && dij < max_rij {
                        bounds.insert((node_i, node_j), dij);
                        bounds.insert((node_j, node_i), dij);
                    } else {
                        bounds.insert((node_i, node_j), bound[1]);
                        bounds.insert((node_j, node_i), bound[1] + dij);
                    }
                } else {
                    if dij > bound[1] && dij < max_rij {
                        bounds.insert((node_i, node_j), dij);
                        bounds.insert((node_j, node_i), max_rij);
                    } else {
                        bounds.insert((node_i, node_j), bound[1]);
                        bounds.insert((node_j, node_i), max_rij);
                    }
                }
            } else {
                bounds.insert((node_i, node_j), bound[1]);
                bounds.insert((node_j, node_i), max_rij);
            }
        }
    }

    bounds
}

/// find some pair of atoms belonging to rigid groups
fn find_rigid_pairs(mol: &Molecule, bounds: &mut Bounds) {
    for a in mol.atoms() {
        for b in mol.atoms() {
            if a.index < b.index {
                let dij = a.distance(b);
                let lij = bounds[&(a.index, b.index)];
                let uij = bounds[&(b.index, a.index)];
                if let Some(nb) = mol.nbonds_between(a.index, b.index) {
                    if nb == 1 {
                        if dij >= lij && dij < uij {
                            bounds.insert((a.index, b.index), dij);
                            bounds.insert((b.index, a.index), dij);
                        } else {
                            let dij = 0.5*(lij + uij);
                            bounds.insert((a.index, b.index), dij);
                            bounds.insert((b.index, a.index), dij);
                        }
                    } else if nb == 2 {
                        if dij >= lij && dij < uij {
                            bounds.insert((a.index, b.index), dij);
                            bounds.insert((b.index, a.index), dij);
                        } else {
                            bounds.insert((a.index, b.index), lij);
                            bounds.insert((b.index, a.index), lij);
                        }
                    }
                }
            }
        }
    }
}

fn get_distance_bounds_v2(mol: &Molecule) -> Bounds {
    let mut bounds = HashMap::new();
    for a in mol.atoms() {
        for b in mol.atoms() {
            if a.index < b.index {
                let dij = a.distance(b);
                let bound = mol.distance_bound(a.index, b.index).unwrap();
                bounds.insert((a.index, b.index), bound[0]);
                bounds.insert((b.index, a.index), bound[1]);
            }
        }
    }

    bounds
}
// 82294367-1b69-4638-a70b-fd8daf02ff3e ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::2351f71f-246f-4193-85c9-7bbe4a9d7587][2351f71f-246f-4193-85c9-7bbe4a9d7587]]
// force component between two atoms
fn get_force_between(lij: f64, uij: f64, dij: f64) -> (f64, f64) {
    let mut force = (lij - dij)/dij;
    let mut weight = 1.0;
    if dij >= lij && dij < uij {
        weight = 0.0;
    } else if dij < lij {
        weight = 1.0;
    } else {
        weight = 1.0;
    }

    (force, weight)
}

impl Molecule {
    pub fn set_momentum(&mut self, index: AtomIndex, m: Point3D) {
        let mut atom = &mut self.graph[index];
        atom.set_momentum(m);
    }

    /// Clean up molecule geometry using stress majorization algorithm
    pub fn clean(&mut self) -> Result<()> {
        let bounds = get_distance_bounds_v1(&self);
        let node_indices: Vec<_> = self.graph.node_indices().collect();
        let nnodes = node_indices.len();

        let maxcycle = nnodes*100;
        let mut icycle = 0;
        let ecut = 1.0E-4;

        let mut ofxz = 0.0;
        loop {
            let mut fxz = 0.0;
            for i in 0..nnodes {
                let node_i = node_indices[i];
                let mut pi = self.get_atom(node_i).expect("atom i from node_i").position();
                let mut disp = [0.0; 3];
                let mut wijs = vec![];
                let npairs = (nnodes - 1) as f64;
                let mut fxzi = 0.0;
                for j in 0..nnodes {
                    if i == j {continue};
                    let node_j = node_indices[j];
                    let pj = self.get_atom(node_j).expect("atom j from node_j").position();
                    let dij = euclidean_distance(pi, pj);
                    let mut bound = [
                        bounds[&(node_i, node_j)],
                        bounds[&(node_j, node_i)],
                    ];
                    if bound[0] > bound[1] {
                        bound.swap(0, 1);
                    }
                    let lij = bound[0];
                    let uij = bound[1];

                    let xij = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
                    let (eij, wij) = get_force_between(lij, uij, dij);
                    let wij = 0.5*wij;
                    wijs.push(wij);
                    let mut fij = [0.0; 3];
                    for v in 0..3 {
                        fij[v] = (1.0 - lij/dij)*xij[v];
                        disp[v] += eij*(pi[v] - pj[v])*wij;
                    }
                    fxzi += wij*(fij[0]*fij[0] + fij[1]*fij[1] + fij[2]*fij[2]);
                }

                let mut swij = 0.0;
                for wij in wijs.iter() {
                    swij += wij;
                }

                // skip updating node_i if all pair weights are zero
                if swij.abs() >= 1e-4 {
                    for v in 0..3 {
                        pi[v] += disp[v]/swij;
                    }
                    self.set_position(node_i, pi);
                }
                // println!("{:?}", (i, fxzi, swij));
                fxz += fxzi;
            }

            println!("cycle: {} energy = {:?}", icycle, fxz);

            if fxz.is_nan() {
                bail!("found invalid number: {:?}", fxz);
            }

            if fxz < ecut || (fxz - ofxz).abs() < ecut || (fxz - ofxz).abs() / fxz < ecut{
                break;
            }

            icycle += 1;
            if icycle > maxcycle {
                break;
            }

            ofxz = fxz;
        }

        Ok(())
    }

}
// 2351f71f-246f-4193-85c9-7bbe4a9d7587 ends here

// fragmentation
// #+name: f0258648-03f4-41c9-949e-f3677c3b44bc

impl Molecule {
    /// Break molecule into smaller fragments based on its bonding connectivity
    pub fn fragment(&self) -> Vec<Molecule> {
        let graph = &self.graph;

        let mut mols = vec![];
        let subgraphs = connected_component_subgraphs(graph);

        for g in subgraphs {
            mols.push(Molecule::from_graph(g));
        }

        mols
    }
}

/// Generate connected components as subgraphs
fn connected_component_subgraphs(graph: &MolGraph) -> Vec<MolGraph>{
    // get fragments from connected components
    let components = petgraph::algo::kosaraju_scc(graph);

    let mut graphs = vec![];
    for nodes in components {
        let g: MolGraph = graph.filter_map
            (
                // node closure:
                // keep atoms in the same component
                |i, n| {
                    if nodes.contains(&i) {
                        Some(n.clone())
                    } else {
                        None
                    }
                },
                // edge closure:
                // keep the bond if bonded atoms are both in the same component
                |i, e| {
                    let (n1, n2) = graph.edge_endpoints(i).unwrap();
                    if nodes.contains(&n1) && nodes.contains(&n2) {
                        Some(e.clone())
                    } else {
                        None
                    }
                }
            );
        graphs.push(g);
    }

    graphs
}

// [[file:~/Workspace/Programming/gchemol/gchemol.note::9b65e04a-b7a8-4118-a2df-d0345c13832c][9b65e04a-b7a8-4118-a2df-d0345c13832c]]
impl Molecule {
    // FIXME: add bonds
    /// Create molecule from small fragments (molecules)
    pub fn combined(mols: &Vec<Molecule>) -> Self {
        let mut mol = Molecule::new("combined");
        for m in mols {
            // add atoms
            for a in m.atoms() {
                let n = mol.add_atom(a.clone());
            }
            // add bonds
            for b in m.bonds() {
                //
            }
        }

        mol
    }
}
// 9b65e04a-b7a8-4118-a2df-d0345c13832c ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::ddf54b1b-6bda-496a-8444-b9762645cc94][ddf54b1b-6bda-496a-8444-b9762645cc94]]
use std::iter::IntoIterator;

fn get_reduced_symbols<'a, I>(symbols: I) -> HashMap<String, usize>
where
    I: IntoIterator,
    I::Item: fmt::Display,
{
    let symbols: Vec<_> = symbols.into_iter()
        .map(|item| format!("{:}", item))
        .collect();

    // 1. count symbols: CCCC ==> C 4
    let mut counts = HashMap::new();
    for x in symbols {
        let c = counts.entry(x).or_insert(0);
        *c += 1;
    }

    counts
}

pub fn get_reduced_formula<'a, I>(symbols: I) -> String
where
    I: IntoIterator,
    I::Item: fmt::Display,
{
    let counts = get_reduced_symbols(symbols);

    let mut syms: Vec<String> = Vec::new();
    let mut to_append = String::new();
    // 2. format the formula
    for (k, v) in counts {
        // 2.1 omit number if the count is 1: C1H4 ==> CH4
        let mut s = String::new();
        if v > 1 {
            s = v.to_string();
        }
        // 2.2 special treatments for C and H
        let reduced = format!("{}{}", k, s);
        if k == "C" {
            syms.insert(0, reduced);
        } else if k == "H" {
            to_append = reduced;
        } else {
            syms.push(reduced);
        }
    }
    // 3. final output
    syms.push(to_append);
    let formula = syms.join("");
    formula
}

impl Molecule {
    /// Return the molecule formula represented in string
    /// Return empty string if molecule containing no atom
    pub fn formula(&self) -> String
    {
        get_reduced_formula(self.symbols())
    }

    // FIXME: add test
    pub fn reduced_symbols(&self) -> HashMap<String, usize> {
        get_reduced_symbols(self.symbols())
    }
}

#[test]
fn test_formula() {
    let symbols   = vec!["C", "H", "C", "H", "H", "H"];
    let formula = get_reduced_formula(&symbols);
    assert_eq!("C2H4", formula);
    let symbols   = vec!["C", "H", "C", "H", "H", "O", "H", "O"];
    let formula = get_reduced_formula(&symbols);
    assert_eq!("C2O2H4", formula);
}
// ddf54b1b-6bda-496a-8444-b9762645cc94 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::5052eafc-f1ab-4612-90d7-0924c3bacb16][5052eafc-f1ab-4612-90d7-0924c3bacb16]]
#[cfg(test)]
mod test {
    use super::*;
    use io;

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
        let a = mol.get_atom(0).expect("failed to get atom with index 0");
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

        mol.set_scaled_positions(&fxyzs);
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
        mol.set_positions(positions.to_vec());
        let a = mol.get_atom(0).expect("failed to get atom with index 0");
        assert_eq!(a.position()[0], -0.90203687);

        // loop over fragments
        let frags = mol.fragment();
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

    #[test]
    fn test_molecule_neighbors() {
        let mol = Molecule::from_file("tests/files/mol2/alanine-gv.mol2").unwrap();
        assert_eq!(12, mol.natoms());

        let atoms = mol.view_atoms();
        let a3 = &atoms[3];

        let ns = a3.neighbors(&mol);
        assert_eq!(4, ns.len());

        let bonds = mol.view_bonds();
        let b35 = &bonds[(3, 5)];

        let (p1, p2) = b35.partners(&mol).unwrap();
    }

    #[test]
    fn test_molecule_pbc_distance() {
        let mut mols = io::read("tests/files/cif/MS-MOR.cif")
            .expect("structure from cif file");
        let mut mol = &mut mols[0];
        let d = mol.distance(AtomIndex::new(0), AtomIndex::new(12))
            .expect("distance between 0 and 12");
        assert_relative_eq!(12.6753, d, epsilon=1e-4);

        // remove periodic bound
        mol.unbuild_crystal();
        let d = mol.distance(AtomIndex::new(0), AtomIndex::new(12))
            .expect("distance between 0 and 12");
        assert_relative_eq!(16.203993, d, epsilon=1e-4);
    }
}
// 5052eafc-f1ab-4612-90d7-0924c3bacb16 ends here
