mod dehydrogenate;
mod matrix;
mod ordering;
mod uniquify;

use std::{collections::BTreeMap, hash::Hash, sync::Arc};

use dehydrogenate::generate_dehydrogenated;
use fxhash::FxHashMap;
use matrix::{HydroCarbonMatrixIter, MatrixFeatures, SymmetricBitMatrix, SymmetricTwoBitsMatrix};
use ordering::RowOrderStore;
use uniquify::gather_unique_matrices_with_symmetry;

#[derive(Debug)]
struct FxHashMapChunkIter<K, V> {
    iter: std::collections::hash_map::IntoIter<K, V>,
    chunk_size: usize,
}

impl<K, V> FxHashMapChunkIter<K, V> {
    fn new(map: FxHashMap<K, V>, chunk_size: usize) -> Self {
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

fn collect_feat2mats<const N: usize>(
    iter: impl IntoIterator<Item = SymmetricBitMatrix<N>>,
) -> FxHashMap<MatrixFeatures<N>, Vec<SymmetricBitMatrix<N>>> {
    let mut hash2mats = FxHashMap::default();
    for mat in iter {
        let (mat, hash) = mat.canonicalize();
        hash2mats.entry(hash).or_insert_with(Vec::new).push(mat);
    }
    hash2mats
}

fn gather_feat2mats<const N: usize>(
    digits: usize,
) -> anyhow::Result<FxHashMap<MatrixFeatures<N>, Vec<SymmetricBitMatrix<N>>>> {
    let mut handlers = Vec::new();
    let iters = HydroCarbonMatrixIter::<N>::new_split(digits)?;
    for iter in iters {
        handlers.push(std::thread::spawn(move || collect_feat2mats(iter)));
    }

    let mut feat2mats = FxHashMap::default();
    for handler in handlers {
        let f2m = handler.join().unwrap();
        for (hash, mut mats) in f2m {
            feat2mats
                .entry(hash)
                .or_insert_with(Vec::new)
                .append(&mut mats);
        }
    }

    Ok(feat2mats)
}

fn generate_hydrocarbons_from_feat2mats<const N: usize>(
    hash2mats: FxHashMap<MatrixFeatures<N>, Vec<SymmetricBitMatrix<N>>>,
    store: &RowOrderStore<N>,
) -> Vec<SymmetricTwoBitsMatrix<N>> {
    let mut all_mats = Vec::new();
    for (hash, mats) in hash2mats {
        let unique_mat_syms = if N <= 11 {
            gather_unique_matrices_with_symmetry::<N, u64>(&mats, &hash, &store)
        } else {
            gather_unique_matrices_with_symmetry::<N, u128>(&mats, &hash, &store)
        };
        for (mat, symmetry) in unique_mat_syms {
            let dehydrogenated = generate_dehydrogenated(mat.into(), &symmetry);
            all_mats.extend(dehydrogenated);
        }
    }
    all_mats
}

pub fn generate_hydrocarbons<const N: usize>(
    digits: usize,
    num_threads: usize,
) -> anyhow::Result<Vec<SymmetricTwoBitsMatrix<N>>> {
    match N {
        0 => return Err(anyhow::anyhow!("N must be greater than 0")),
        17.. => return Err(anyhow::anyhow!("N must be less than or equal to 16")),
        _ => (),
    };

    let feat2mats = gather_feat2mats::<N>(digits)?;
    let store = Arc::new(RowOrderStore::<N>::new_parallel(num_threads as u64));

    let mut handlers = Vec::new();
    let chunk_size = feat2mats.len().div_ceil(num_threads);
    for sub_feat2mats in FxHashMapChunkIter::new(feat2mats, chunk_size) {
        let store = store.clone();
        handlers.push(std::thread::spawn(move || {
            generate_hydrocarbons_from_feat2mats(sub_feat2mats, &store)
        }));
    }

    let mats = handlers
        .into_iter()
        .flat_map(|h| h.join().unwrap())
        .collect();
    Ok(mats)
}

pub fn display_counts<const N: usize>(mats: &[SymmetricTwoBitsMatrix<N>]) {
    let mut num_h_to_count = BTreeMap::new();
    for mat in mats {
        let num_h = mat.count_hydrogen();
        *num_h_to_count.entry(num_h).or_insert(0) += 1;
    }

    println!("===== [C = {N:>2}] =====");
    println!("#H: #HydroCarbods");
    for (num_h, count) in num_h_to_count {
        println!("{:>2}: {}", num_h, count);
    }
}
