use std::{io, path::PathBuf, sync::Arc};

use dehydrogenate::generate_all_dehydrogenated;
use fxhash::FxHashMap;
use itertools::Itertools;
use matrix::{get_symmetry, load_matrices};
use search::RowOrderStore;

mod dehydrogenate;
mod matrix;
mod search;

const NUM_THREADS: usize = 128;

fn do_generate<const N: usize>() -> io::Result<()> {
    let root = PathBuf::from("input");
    let path = root.join(format!("C{N:0>2}.json"));
    let base_matrices = load_matrices::<N>(&path)?;
    let store = Arc::new(RowOrderStore::<N>::new());

    let mut handlers = Vec::new();
    let chunk_size = base_matrices.len().div_ceil(NUM_THREADS);
    for sub_matrices in base_matrices
        .into_iter()
        .chunks(chunk_size)
        .into_iter()
        .map(|chunk| chunk.collect_vec())
    {
        let store = store.clone();
        let handler = std::thread::spawn(move || {
            let mut sub_mats = Vec::new();
            for (mat, hash) in sub_matrices {
                let symmetry = get_symmetry(&mat, &hash, &store);
                sub_mats.extend(generate_all_dehydrogenated(mat, &symmetry));
            }
            sub_mats
        });
        handlers.push(handler);
    }
    let all_mats = handlers
        .into_iter()
        .flat_map(|handler| handler.join().unwrap())
        .collect_vec();

    let mut num_h_to_mats = FxHashMap::default();
    for mat in all_mats {
        let num_h = mat.count_hydrogen();
        num_h_to_mats
            .entry(num_h)
            .or_insert_with(Vec::new)
            .push(mat);
    }
    let mut keys = num_h_to_mats.keys().copied().collect_vec();
    keys.sort();

    println!("===== [C = {N}] =====");
    println!("#H: # of structures");
    for key in keys {
        println!("{:>2}: {}", key, num_h_to_mats[&key].len());
    }
    Ok(())
}

macro_rules! do_generate_many {
    ($($n:expr),+) => {
        $(
            do_generate::<$n>()?;
        )*
    };
}

fn main() -> io::Result<()> {
    do_generate_many!(2, 3, 4, 5, 6, 7, 8, 9, 10);
    Ok(())
}
