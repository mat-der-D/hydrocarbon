use std::{collections::VecDeque, hash::Hash};

use fxhash::{FxHashMap, FxHashSet};

use super::{
    matrix::{MatrixFeatures, SymmetricBitMatrix},
    ordering::{Permutation, PermutationsStore},
};

pub fn gather_unique_matrices_with_symmetry<const N: usize, T>(
    matrices: &[SymmetricBitMatrix<N>],
    feature: &MatrixFeatures<N>,
    store: &PermutationsStore<N>,
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
    store: &PermutationsStore<N>,
) -> (FxHashSet<T>, Vec<Permutation<N>>)
where
    T: From<SymmetricBitMatrix<N>> + Copy + Eq + Hash,
{
    let mut symmetry = FxHashSet::default();

    let mut hash2perm = FxHashMap::default();
    hash2perm.insert(hash, Permutation::<N>::IDENTITY);

    let mut queue = VecDeque::new();
    queue.push_back((*mat, Permutation::<N>::IDENTITY));

    let generators = feature.generate_permutations(store);
    while let Some((mat_, perm)) = queue.pop_front() {
        for gen in generators {
            let new_mat = gen.permute(&mat_);
            let new_mat_hash: T = new_mat.into();
            let new_mat_perm = *gen * perm;
            if let Some(saved_perm) = hash2perm.get(&new_mat_hash) {
                if new_mat_perm != *saved_perm {
                    symmetry.insert(saved_perm.inverse() * new_mat_perm);
                }
            } else {
                hash2perm.insert(new_mat_hash, new_mat_perm);
                queue.push_back((new_mat, new_mat_perm));
            }
        }
    }
    let variants = hash2perm.into_keys().collect();
    (variants, symmetry.into_iter().collect())
}
