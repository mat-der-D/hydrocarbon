use std::{collections::BTreeMap, sync};

use dehydrogenate::generate_all_dehydrogenated;
use fxhash::FxHashMap;
use itertools::Itertools;
use search::{make_unique, MatrixSearcher, RowOrderStore};

mod dehydrogenate;
mod matrix;
mod search;

fn generate_hydrocarbons<const N: usize>(num_threads: usize) {
    if N > 16 {
        panic!("N must be less than or equal to 16");
    }
    let mut hash2mat = FxHashMap::default();
    let mut searcher = MatrixSearcher::<N>::new();
    searcher.search(|mat| {
        let mat = mat.partially_canonicalize();
        let hash = mat.make_hash();
        hash2mat.entry(hash).or_insert_with(Vec::new).push(mat);
    });

    let store = sync::Arc::new(RowOrderStore::<N>::new());
    let mut handlers = Vec::new();
    let chunk_size = hash2mat.len().div_ceil(num_threads);
    for sub_hash2mat in hash2mat
        .into_iter()
        .chunks(chunk_size)
        .into_iter()
        .map(|c| c.collect::<FxHashMap<_, _>>())
    {
        let store = store.clone();
        let handler = std::thread::spawn(move || {
            let mut all_mats = Vec::new();
            for (hash, mats) in sub_hash2mat {
                let unique_mat_syms = make_unique(&mats, &hash, &store);
                for (mat, symmetry) in unique_mat_syms {
                    let dehydrogenated = generate_all_dehydrogenated(mat.into(), &symmetry);
                    all_mats.extend(dehydrogenated);
                }
            }
            all_mats
        });
        handlers.push(handler);
    }
    let all_mats = handlers
        .into_iter()
        .flat_map(|h| h.join().unwrap())
        .collect_vec();

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
    generate_hydrocarbons_many!(2, 3, 4, 5, 6, 7, 8, 9, 10);
}
