mod matrix;
mod multi_bond;
mod ordering;
mod single_bond;

use std::{collections::BTreeMap, hash::Hash, sync::Arc};

use fxhash::FxHashMap;
use matrix::{HydroCarbonMatrixIter, MatrixHash, SymmetricBitMatrix};
use multi_bond::generate_all_dehydrogenated;
use ordering::RowOrderStore;
use single_bond::make_unique;

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

fn gather_hash2mats_single_thread<const N: usize>(
) -> FxHashMap<MatrixHash<N>, Vec<SymmetricBitMatrix<N>>> {
    let iter = HydroCarbonMatrixIter::<N>::default();
    collect_hash2mats(iter)
}

fn collect_hash2mats<const N: usize>(
    iter: impl IntoIterator<Item = SymmetricBitMatrix<N>>,
) -> FxHashMap<MatrixHash<N>, Vec<SymmetricBitMatrix<N>>> {
    let mut hash2mats = FxHashMap::default();
    for mat in iter {
        let (mat, hash) = mat.canonicalize();
        hash2mats.entry(hash).or_insert_with(Vec::new).push(mat);
    }
    hash2mats
}

fn gather_hash2mats<const N: usize>(
    digits: usize,
) -> FxHashMap<MatrixHash<N>, Vec<SymmetricBitMatrix<N>>> {
    let mut handlers = Vec::new();
    let iters = HydroCarbonMatrixIter::<N>::create_multi(digits);
    for iter in iters {
        handlers.push(std::thread::spawn(move || collect_hash2mats(iter)));
    }

    let mut hash2mats = FxHashMap::default();
    for handler in handlers {
        let h2m = handler.join().unwrap();
        for (hash, mut mats) in h2m {
            hash2mats
                .entry(hash)
                .or_insert_with(Vec::new)
                .append(&mut mats);
        }
    }

    hash2mats
}

fn generate_hydrocarbons<const N: usize>(num_threads: usize, num_digits: usize) {
    let hash2mat = match N {
        0 => panic!("N must be greater than 0"),
        17.. => panic!("N must be less than or equal to 16"),
        1..=3 => gather_hash2mats_single_thread::<N>(),
        _ => gather_hash2mats::<N>(num_digits),
    };
    let store = Arc::new(RowOrderStore::<N>::new_parallel(num_threads as u64));

    let mut handlers = Vec::new();
    let chunk_size = hash2mat.len().div_ceil(num_threads);
    for sub_hash2mat in FxHashMapChunkIter::new(hash2mat, chunk_size) {
        let store = store.clone();
        let handler = std::thread::spawn(move || {
            let mut all_mats = Vec::new();
            for (hash, mats) in sub_hash2mat {
                let unique_mat_syms = if N <= 11 {
                    make_unique::<N, u64>(&mats, &hash, &store)
                } else {
                    make_unique::<N, u128>(&mats, &hash, &store)
                };
                for (mat, symmetry) in unique_mat_syms {
                    let dehydrogenated = generate_all_dehydrogenated(mat.into(), &symmetry);
                    all_mats.extend(dehydrogenated);
                }
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
    ($($n:expr),+;$d:expr) => {
        $(
            generate_hydrocarbons::<$n>(NUM_THREADS, $d);
        )*
    };
}

fn main() {
    let digits = 7;
    generate_hydrocarbons_many!(2, 3, 4, 5, 6, 7, 8, 9, 10; digits); // TODO: N = 1 に対応する
}
