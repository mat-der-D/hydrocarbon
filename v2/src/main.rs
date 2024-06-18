use fxhash::FxHashMap;
use matrix::HydroCarbonMatrixIter;

mod matrix;

fn main() {
    const N: usize = 10;
    let instant = std::time::Instant::now();
    let mut hash2mats = FxHashMap::default();
    for mat in HydroCarbonMatrixIter::<N>::new() {
        let (mat, hash) = mat.canonicalize();
        hash2mats.entry(hash).or_insert_with(Vec::new).push(mat);
    }
    println!("# of keys: {}", hash2mats.len());
    println!("{:?}", instant.elapsed());
}
