// [[file:~/Workspace/Programming/gchemol/core/gchemol-core.note::45c615d8-4272-4f7a-a196-71f3cc1c4bef][45c615d8-4272-4f7a-a196-71f3cc1c4bef]]
/// search neighbors for aperiodic system
pub fn neighbors_for_aperiodic
    (
        positions: &Vec<[f64; 3]>,
        cutoff: f64
    )
    -> HashMap<usize, Vec<(usize, f64, Vector3<f64>)>>
{
    let tree = Octree::new(&positions);

    let mut kneighbors = HashMap::new();
    // run queries over all relevant images
    for (current, &p) in positions.iter().enumerate() {
        let mut neighbors = Vec::new();
        let nps = tree.search(p, cutoff);
        for &(index, distance) in nps.iter() {
            if index != current {
                neighbors.push((index, distance, Vector3::new(0.0, 0.0, 0.0)));
            }
        }
        kneighbors.insert(current, neighbors);
    }

    kneighbors
}

/// search neighbors for periodic system
pub fn neighbors_for_periodic
    (
        positions: &Vec<[f64; 3]>,
        cell: UnitCell, cutoff: f64
    )
    -> HashMap<usize, Vec<(usize, f64, Vector3<f64>)>>
{
    let images = cell.relevant_images(cutoff);

    let tree = Octree::new(&positions);

    let mut kneighbors = HashMap::new();
    // run queries over all relevant images
    for (current, &p) in positions.iter().enumerate() {
        // to avoid octree building for each image
        // we mirror the query points and then mirror back
        let mut neighbors: Vec<(usize, f64, Vector3<f64>)> = Vec::new();
        for &image in images.iter() {
            let disp = cell.matrix * image;
            let vquery = Vector3::from(p) + disp;
            let nps = tree.search(*vquery.as_ref(), cutoff);
            for &(index, distance) in nps.iter() {
                if index != current {
                    neighbors.push((index, distance, -1.*image))
                }
            }
        }
        kneighbors.insert(current, neighbors);
    }

    kneighbors
}
// 45c615d8-4272-4f7a-a196-71f3cc1c4bef ends here
