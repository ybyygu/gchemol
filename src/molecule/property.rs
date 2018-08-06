// [[file:~/Workspace/Programming/gchemol/gchemol.note::a2f1dbf4-c41d-4ca4-986e-5e8cfb3d5d08][a2f1dbf4-c41d-4ca4-986e-5e8cfb3d5d08]]
use std::collections::HashMap;
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
    pub fn new() -> Self {
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
