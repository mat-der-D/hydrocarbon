use std::{fs, io, path, sync};

use fxhash::FxHashMap;
use itertools::Itertools;

mod matrix;
mod search;

use matrix::{MatrixAndHash, SymmetricBitMatrix};
use search::{make_unique, MatrixSearcher, RowOrderStore};

fn do_search<const N: usize>(num_threads: usize, save_results: bool) -> io::Result<()> {
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
            let mut unique_matrices = Vec::new();
            for (hash, mats) in sub_hash2mat.iter() {
                unique_matrices.extend(make_unique(mats, hash, &store));
            }
            unique_matrices
        });
        handlers.push(handler);
    }

    let unique_mats = handlers
        .into_iter()
        .flat_map(|h| h.join().unwrap())
        .collect::<Vec<_>>();

    println!("{:>2} ==> {:>6}", N, unique_mats.len());
    if !save_results {
        return Ok(());
    }

    let root_dir = path::PathBuf::from("output");
    fs::create_dir_all(&root_dir)?;
    let file_path = root_dir.join(format!("C{:0>2}.json", N));
    let mut writer = io::BufWriter::new(fs::File::create(file_path)?);
    write_matrices(&unique_mats, &mut writer)?;
    Ok(())
}

fn write_matrices<const N: usize>(
    matrices: &[SymmetricBitMatrix<N>],
    writer: &mut impl io::Write,
) -> io::Result<()> {
    let mat_hash_vec = matrices
        .iter()
        .map(|mat| MatrixAndHash::new(*mat))
        .collect_vec();
    serde_json::to_writer_pretty(writer, &mat_hash_vec)?;
    Ok(())
}

const SAVE_RESULTS: bool = false; // 結果を保存する場合はここを true にする
const NUM_THREADS: usize = 128;

macro_rules! do_search_many {
    ($($mat_size:expr),+) => {
        $(
            do_search::<$mat_size>(NUM_THREADS, SAVE_RESULTS)?;
        )*
    };
}

fn main() -> io::Result<()> {
    println!("===== # of hydrocarbons with only single-bonds =====");
    println!("#C ==> # of structures");
    do_search_many!(2, 3, 4, 5, 6, 7, 8, 9, 10);
    Ok(())
}
