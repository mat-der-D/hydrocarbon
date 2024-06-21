use fxhash::FxHashSet;

use super::matrix::SymmetricTwoBitsMatrix;

pub fn generate_dehydrogenated<const N: usize>(
    base: SymmetricTwoBitsMatrix<N>,
    symmetry: &[[usize; N]],
) -> Vec<SymmetricTwoBitsMatrix<N>> {
    let mut mats = vec![base];
    let base_set: FxHashSet<_> = [base].into_iter().collect();
    let mut seen = vec![base_set];
    let mut current = vec![base];
    let mut next = Vec::new();

    while !current.is_empty() {
        for mat in current {
            for (row, col) in mat.dehydrogenatable_bonds() {
                let new_mat = mat.create_dehydrogenated_unchecked(row, col);
                if seen.iter().any(|s| s.contains(&new_mat)) {
                    continue;
                }
                mats.push(new_mat);
                next.push(new_mat);

                let variants = make_variants(&new_mat, symmetry);
                seen.push(variants);
            }
        }

        current = next;
        next = Vec::new();
        seen.clear();
    }
    mats
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
