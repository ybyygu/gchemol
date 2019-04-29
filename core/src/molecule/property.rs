// property.rs
// :PROPERTIES:
// :header-args: :tangle src/molecule/property.rs
// :END:

// [[file:~/Workspace/Programming/gchemol-rs/gchemol/core/gchemol-core.note::*property.rs][property.rs:1]]
use std::collections::HashMap;

use serde::{Deserialize, Serialize, de::DeserializeOwned};
// use serde::de::DeserializeOwned;
//use serde::{de::DeserializeOwned, ser::Serialize};

use serde_json;
use std::result;

/// A container storing extra information managed as key/value pairs
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PropertyStore {
    data: HashMap<String, String>,
}

impl PropertyStore {
    pub fn new() -> Self {
        PropertyStore {
            data: HashMap::new(),
        }
    }

    /// Returns true if the map contains a value for the specified key.
    pub fn contains_key(&self, key: &str) -> bool {
        self.data.contains_key(key)
    }

    /// Retrieve property associated with the `key`
    pub fn load<D: DeserializeOwned>(&self, key: &str) -> result::Result<D, serde_json::Error> {
        let serialized = self.data.get(key).unwrap();
        serde_json::from_str(&serialized)
    }

    /// Store property associatd with the `key`
    pub fn store<D: Serialize>(&mut self, key: &str, value: D) {
        let serialized = serde_json::to_string(&value).unwrap();
        self.data.insert(key.into(), serialized);
    }

    /// Discard property associated with the `key`
    pub fn discard(&mut self, key: &str) {
        self.data.remove(key);
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
// property.rs:1 ends here
