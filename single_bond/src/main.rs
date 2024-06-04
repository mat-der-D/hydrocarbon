use std::sync::Arc;

use fxhash::FxHashMap;
use itertools::Itertools;

mod matrix;
mod search;

use search::{make_unique, MatrixSearcher, RowOrderStore};

const MAT_SIZE: usize = 10;
const NUM_THREADS: usize = 128;

fn main() {
    let mut hash2mat = FxHashMap::default();
    let mut searcher = MatrixSearcher::<MAT_SIZE>::new();
    searcher.search(|mat| {
        let mat = mat.partially_canonicalize();
        let hash = mat.make_hash();
        hash2mat.entry(hash).or_insert_with(Vec::new).push(mat);
    });
    println!("# of Keys: {}", hash2mat.len());

    let store = Arc::new(RowOrderStore::<MAT_SIZE>::new());
    let mut handlers = Vec::new();
    let chunk_size = hash2mat.len().div_ceil(NUM_THREADS);
    for (i, sub_hash2mat) in hash2mat
        .into_iter()
        .chunks(chunk_size)
        .into_iter()
        .map(|c| c.collect::<FxHashMap<_, _>>())
        .enumerate()
    {
        let store = store.clone();
        let handler = std::thread::spawn(move || {
            println!("Thread {:0>3} started", i);
            let mut unique_matrices = Vec::new();
            for (hash, mats) in sub_hash2mat.iter() {
                unique_matrices.extend(make_unique(mats, hash, &store));
            }
            println!("Thread {:0>3} finished", i);
            unique_matrices
        });
        handlers.push(handler);
    }

    let unique_mats = handlers
        .into_iter()
        .flat_map(|h| h.join().unwrap())
        .collect::<Vec<_>>();

    println!("# of Unique matrices: {}", unique_mats.len());
}
