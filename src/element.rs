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
#[derive(Debug, Clone)]
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
            &Element(num) => ELEMENT_DATA[num-1].1.to_string(),
            &Dummy(ref sym) => format!("dummy atom {}", sym),
        }
    }
}

impl fmt::Display for AtomKind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            &Element(num) => write!(f, "Element {:}", self.symbol()),
            &Dummy(ref sym) =>  write!(f, "Dummy atom {:}", self.symbol()),
        }
    }
}
// 731651b9-ebba-4df3-86f7-083f837e4065 ends here

// [[file:~/Workspace/Programming/gchemol/gchemol.note::b95edc21-e696-4625-ba99-94257394772d][b95edc21-e696-4625-ba99-94257394772d]]
use self::AtomKind::{Element, Dummy};

/// Return AtomKind using common sense
pub fn atom_kind_from_string<T: Into<String>>(sym: T) -> AtomKind {
    let sym = sym.into();
    for (i, &(s, n) ) in ELEMENT_DATA.iter().enumerate() {
        if s == sym  || n == sym {
            return Element(i+1);
        }
    }

    Dummy(sym)
}

#[test]
fn test_element() {
    let x = Element(12);
    println!("{:}", x);
    assert_eq!(12, x.number());
    assert_eq!("Mg", x.symbol());
    assert_eq!("magnesium", x.name());

    let x = Dummy("X".to_string());
    assert_eq!("X", x.symbol());
    assert_eq!(0, x.number());
}
// b95edc21-e696-4625-ba99-94257394772d ends here
