// [[file:~/Workspace/Programming/gchemol/gchemol.note::618eea82-a702-4d2f-a873-3807ead50d4b][618eea82-a702-4d2f-a873-3807ead50d4b]]
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::error::Error;

use Point3D;
use Points;
/// write coordinates in xyz format
pub fn write_as_xyz(points: &Points, filename: &str) -> Result<(), Box<Error>>
{
    let mut lines = String::new();

    // meta information
    lines.push_str(format!("{}\n", points.len()).as_str());
    lines.push_str("default title\n");

    // coordinates
    for &p in points.iter() {
        let s = format!("C {:-18.6}{:-18.6}{:-18.6}\n", p[0], p[1], p[2]);
        lines.push_str(s.as_str());
    }

    // final empty line
    lines.push_str("\n");

    // save as a file
    let f = File::create(filename)?;
    let mut writer = BufWriter::new(&f);
    writer.write_all(&lines.as_bytes())?;

    Ok(())
}
// 618eea82-a702-4d2f-a873-3807ead50d4b ends here
