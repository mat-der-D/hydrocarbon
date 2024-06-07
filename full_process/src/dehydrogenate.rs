use fxhash::FxHashSet;

use crate::matrix::SymmetricTwoBitsMatrix;

pub fn generate_all_dehydrogenated<const N: usize>(
    base_matrix: SymmetricTwoBitsMatrix<N>,
    symmetry: &[[usize; N]],
) -> Vec<SymmetricTwoBitsMatrix<N>> {
    let mut all_mats = vec![base_matrix];
    let mut seen = vec![[base_matrix].iter().copied().collect::<FxHashSet<_>>()];
    let mut current = vec![base_matrix];
    let mut next = Vec::new();

    while !current.is_empty() {
        for mat in current {
            for (row, col) in mat.dehydrogenatable_bonds() {
                let mat_new = mat.create_dehydrogenated_unchecked(row, col);
                if seen
                    .iter()
                    .any(|s: &FxHashSet<SymmetricTwoBitsMatrix<N>>| s.contains(&mat_new))
                {
                    continue;
                }
                all_mats.push(mat_new);
                next.push(mat_new);

                let variants = make_variants(&mat_new, symmetry);
                seen.push(variants);
            }
        }

        current = next;
        next = Vec::new();
        seen.clear();
    }
    all_mats
}

fn make_variants<const N: usize>(
    mat: &SymmetricTwoBitsMatrix<N>,
    symmetry: &[[usize; N]],
) -> FxHashSet<SymmetricTwoBitsMatrix<N>> {
    symmetry
        .iter()
        .map(|row_order| mat.create_rearranged(row_order))
        .collect()
}
