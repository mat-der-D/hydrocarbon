use fxhash::FxHashSet;

use crate::matrix::{MatrixHash, SymmetricTwoBitsMatrix};

pub fn generate_all_dehydrogenated2<const N: usize>(
    base_matrix: SymmetricTwoBitsMatrix<N>,
    hash: &MatrixHash<N>,
) -> Vec<SymmetricTwoBitsMatrix<N>> {
    let mut all_mats = vec![base_matrix];
    let mut seen = vec![[base_matrix].iter().copied().collect::<FxHashSet<_>>()];
    let mut current = vec![base_matrix];
    let mut next = Vec::new();

    while !current.is_empty() {
        for mat in current {
            for (row, col) in mat.dehydrogenatable_bonds() {
                let mat_new = mat.create_dehydrogenated_unchecked(row, col);
                if all_mats.iter().any(|m| m.is_equivalent_to(&mat_new, hash)) {
                    continue;
                }
                all_mats.push(mat_new);
                next.push(mat_new);
            }
        }

        current = next;
        next = Vec::new();
        seen.clear();
    }
    all_mats
}
