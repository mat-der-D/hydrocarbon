use std::hash::Hash;

use fxhash::FxHashSet;

use crate::{
    matrix::{MatrixHash, SymmetricBitMatrix},
    ordering::RowOrderStore,
};

pub fn make_unique<const N: usize, T>(
    matrices: &[SymmetricBitMatrix<N>],
    hash: &MatrixHash<N>,
    store: &RowOrderStore<N>,
) -> Vec<(SymmetricBitMatrix<N>, Vec<[usize; N]>)>
where
    T: From<SymmetricBitMatrix<N>> + Eq + Hash,
{
    let mut unique_matrix_symmetry_pairs = Vec::new();
    let mut seen = Vec::new();
    for &mat in matrices {
        let mat_hash: T = mat.into();
        if seen.iter().any(|s: &FxHashSet<T>| s.contains(&mat_hash)) {
            continue;
        }

        let (variants, symmetry) = make_variants_symmetry(&mat, mat_hash, hash, store);
        seen.push(variants);
        unique_matrix_symmetry_pairs.push((mat, symmetry));
    }
    unique_matrix_symmetry_pairs
}

fn make_variants_symmetry<const N: usize, T>(
    &mat: &SymmetricBitMatrix<N>,
    mat_hash: T,
    hash: &MatrixHash<N>,
    store: &RowOrderStore<N>,
) -> (FxHashSet<T>, Vec<[usize; N]>)
where
    T: From<SymmetricBitMatrix<N>> + Eq + Hash,
{
    let mut variants = FxHashSet::default();
    let mut symmetry = Vec::new();
    for &row_order in hash.generate_row_orders(store) {
        let rearranged_hash: T = mat.create_rearranged(&row_order).into();
        if mat_hash == rearranged_hash {
            symmetry.push(row_order);
        }
        variants.insert(rearranged_hash);
    }
    (variants, symmetry)
}
