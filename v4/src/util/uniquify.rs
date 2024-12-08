use std::{collections::VecDeque, hash::Hash};

use fxhash::FxHashSet;

use super::{
    matrix::{MatrixFeatures, SymmetricBitMatrix},
    ordering::{Permutation, RowOrderStore},
};

pub fn gather_unique_matrices_with_symmetry<const N: usize, T>(
    matrices: &[SymmetricBitMatrix<N>],
    feature: &MatrixFeatures<N>,
    store: &RowOrderStore<N>,
) -> Vec<(SymmetricBitMatrix<N>, Vec<Permutation<N>>)>
where
    T: From<SymmetricBitMatrix<N>> + Copy + Eq + Hash,
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
) -> (FxHashSet<T>, Vec<Permutation<N>>)
where
    T: From<SymmetricBitMatrix<N>> + Copy + Eq + Hash,
{
    let mut variants = FxHashSet::default();
    variants.insert(hash);
    let mut symmetry = Vec::new();
    let generators = feature.generate_permutations(store);
    let mut queue = VecDeque::new();
    queue.push_back((*mat, Permutation::<N>::identity()));
    while !queue.is_empty() {
        let (mat_, perm) = queue.pop_front().unwrap();
        for gen in generators {
            let new_mat = gen.permute(&mat_);
            let new_mat_hash: T = new_mat.into();
            let accumulated = *gen * perm;
            if hash == new_mat_hash
                && !accumulated.is_identity()
                && !symmetry.contains(&accumulated)
            {
                symmetry.push(accumulated);
            }
            if !variants.contains(&new_mat_hash) {
                variants.insert(new_mat_hash);
                queue.push_back((new_mat, accumulated));
            }
        }
    }
    (variants, symmetry)
}
