use std::hash::Hash;

use fxhash::FxHashSet;

use crate::{
    matrix::{MatrixFeatures, SymmetricBitMatrix},
    ordering::RowOrderStore,
};

pub fn gather_unique_matrices_with_symmetry<const N: usize, T>(
    matrices: &[SymmetricBitMatrix<N>],
    feature: &MatrixFeatures<N>,
    store: &RowOrderStore<N>,
) -> Vec<(SymmetricBitMatrix<N>, Vec<[usize; N]>)>
where
    T: From<SymmetricBitMatrix<N>> + Eq + Hash,
{
    let mut mat_sym_vec = Vec::new();
    let mut seen: Vec<FxHashSet<T>> = Vec::new();
    for &mat in matrices {
        let hash: T = mat.into();
        if seen.iter().any(|s| s.contains(&hash)) {
            continue;
        }

        let (variants, symmetry) = make_variants_symmetry(&mat, hash, feature, store);
        seen.push(variants);
        mat_sym_vec.push((mat, symmetry));
    }
    mat_sym_vec
}

fn make_variants_symmetry<const N: usize, T>(
    mat: &SymmetricBitMatrix<N>,
    hash: T,
    feature: &MatrixFeatures<N>,
    store: &RowOrderStore<N>,
) -> (FxHashSet<T>, Vec<[usize; N]>)
where
    T: From<SymmetricBitMatrix<N>> + Eq + Hash,
{
    let mut variants = FxHashSet::default();
    let mut symmetry = Vec::new();
    for &row_order in feature.generate_row_orders(store) {
        let rearranged_hash: T = mat.create_rearranged(&row_order).into();
        if hash == rearranged_hash {
            symmetry.push(row_order);
        }
        variants.insert(rearranged_hash);
    }
    (variants, symmetry)
}
