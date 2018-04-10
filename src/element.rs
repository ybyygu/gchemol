// [[file:~/Workspace/Programming/gchemol/gchemol.note::731651b9-ebba-4df3-86f7-083f837e4065][731651b9-ebba-4df3-86f7-083f837e4065]]
use std::fmt::{self, Debug, Display};

#[derive(Debug, Copy, Clone)]
pub enum AtomKind<'a> {
    /// physical elements
    Element(usize),
    /// a dummy-atom is not a real atom
    Dummy(&'a str),
}

static ELEMENT_DATA: [(&'static str, &'static str); 118] = [
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


impl <'a>AtomKind<'a> {
    pub fn symbol(&self) -> String {
        match(self) {
            &Element(num) => ELEMENT_DATA[num-1].0.to_string(),
            &Dummy(sym) => sym.to_string()
        }
    }

    pub fn number(&self) -> usize {
        match(self) {
            &Element(num) => num,
            &Dummy(sym) => 0,
        }
    }

    pub fn name(&self) -> String {
        match(self) {
            &Element(num) => ELEMENT_DATA[num-1].1.to_string(),
            &Dummy(sym) => format!("dummy atom {}", sym),
        }
    }
}

impl <'a>fmt::Display for AtomKind<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match(self) {
            &Element(num) => write!(f, "Element {:}", self.symbol()),
            &Dummy(sym) =>  write!(f, "Dummy atom {:}", self.symbol()),
        }
    }
}

use self::AtomKind::{Element, Dummy};

#[test]
fn test_element() {
    let x = Element(12);
    println!("{:}", x);
    println!("symbol = {:}", x.symbol());
    println!("number = {:}", x.number());
    println!("name = {:}", x.name());
    let x = Dummy("X");
    println!("{:}", x);
    println!("symbol = {:}", x.symbol());
    println!("number = {:}", x.number());
    println!("name = {:}", x.name());
}
// 731651b9-ebba-4df3-86f7-083f837e4065 ends here
