use std::{collections::BTreeMap, hash::Hash, sync};

use dehydrogenate::generate_all_dehydrogenated;
use fxhash::FxHashMap;
use search::{make_symmetry, MatrixSearcher, RowOrderStore};

mod dehydrogenate;
mod matrix;
mod search;

#[derive(Debug)]
pub struct FxHashMapChunkIter<K, V> {
    iter: std::collections::hash_map::IntoIter<K, V>,
    chunk_size: usize,
}

impl<K, V> FxHashMapChunkIter<K, V> {
    pub fn new(map: FxHashMap<K, V>, chunk_size: usize) -> Self {
        Self {
            iter: map.into_iter(),
            chunk_size,
        }
    }
}

impl<K, V> Iterator for FxHashMapChunkIter<K, V>
where
    K: Eq + Hash,
{
    type Item = FxHashMap<K, V>;

    fn next(&mut self) -> Option<Self::Item> {
        let map: FxHashMap<K, V> = (&mut self.iter).take(self.chunk_size).collect();
        if map.is_empty() {
            None
        } else {
            Some(map)
        }
    }
}

fn generate_hydrocarbons<const N: usize>(num_threads: usize) {
    if N > 16 {
        panic!("N must be less than or equal to 16");
    }
    let mut hash2mat = FxHashMap::default();
    let mut searcher = MatrixSearcher::<N>::new();
    searcher.search(|mat| {
        let (mat, hash) = mat.partially_canonicalize();
        hash2mat.entry(hash).or_insert(mat);
    });

    let store = sync::Arc::new(RowOrderStore::<N>::new());
    let mut handlers = Vec::new();
    let chunk_size = hash2mat.len().div_ceil(num_threads);
    for sub_hash2mat in FxHashMapChunkIter::new(hash2mat, chunk_size) {
        let store = store.clone();
        let handler = std::thread::spawn(move || {
            let mut all_mats = Vec::new();
            for (hash, mat) in sub_hash2mat {
                let symmetry = make_symmetry(&mat, &hash, &store);
                let dehydrogenated = generate_all_dehydrogenated(mat.into(), &symmetry);
                all_mats.extend(dehydrogenated);
            }
            all_mats
        });
        handlers.push(handler);
    }
    let all_mats: Vec<_> = handlers
        .into_iter()
        .flat_map(|h| h.join().unwrap())
        .collect();

    // 結果表示用の集計
    let mut num_h_to_count = BTreeMap::new();
    for mat in all_mats {
        let num_h = mat.count_hydrogen();
        *num_h_to_count.entry(num_h).or_insert(0u32) += 1;
    }

    // 結果表示
    println!("===== [C = {N:>2}] =====");
    println!("#H: # of structures");
    for (num_h, count) in num_h_to_count {
        println!("{:>2}: {}", num_h, count);
    }
}

const NUM_THREADS: usize = 128;

macro_rules! generate_hydrocarbons_many {
    ($($n:expr),+) => {
        $(
            generate_hydrocarbons::<$n>(NUM_THREADS);
        )*
    };
}

fn main() {
    generate_hydrocarbons_many!(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
}
