use crate::matrix::{MatrixHash, SymmetricBitMatrix};

#[derive(Debug)]
pub struct MatrixSearcher<const N: usize> {
    matrix: SymmetricBitMatrix<N>,
}

impl<const N: usize> Default for MatrixSearcher<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> MatrixSearcher<N> {
    pub fn new() -> Self {
        Self {
            matrix: SymmetricBitMatrix::zero(),
        }
    }

    pub fn search(&mut self, mut f: impl FnMut(&SymmetricBitMatrix<N>)) {
        self.matrix = SymmetricBitMatrix::zero();
        self.search_impl(0, 1, &mut f);
    }

    fn search_impl(&mut self, row: usize, col: usize, f: &mut impl FnMut(&SymmetricBitMatrix<N>)) {
        if row == N - 1 {
            f(&self.matrix);
            return;
        }

        let (next_row, next_col) = if col == N - 1 {
            (row + 1, row + 2)
        } else {
            (row, col + 1)
        };

        if Self::is_fine_up_to(&self.matrix, next_row) {
            self.search_impl(next_row, next_col, f);
        }
        self.matrix.activate(row, col);
        if Self::is_fine_up_to(&self.matrix, next_row) {
            self.search_impl(next_row, next_col, f);
        }
        self.matrix.deactivate(row, col);
    }

    fn is_fine_up_to(matrix: &SymmetricBitMatrix<N>, idx_row: usize) -> bool {
        let mut prev_row = u16::MAX;
        let mut prev_nbits = 4;
        for (i, &row) in matrix.rows().iter().enumerate() {
            let nbits = row.count_ones();
            if nbits > prev_nbits || (nbits == prev_nbits && row > prev_row) {
                return false;
            }
            if i < idx_row {
                prev_row = row;
                prev_nbits = nbits;
            }
        }
        idx_row < N - 1 || matrix.is_connected()
    }
}

pub fn make_unique2<const N: usize>(
    matrices: &[SymmetricBitMatrix<N>],
    hash: &MatrixHash<N>,
) -> Vec<SymmetricBitMatrix<N>> {
    let mut unique_matrices = Vec::new();
    let mut seen = Vec::new();
    for &mat in matrices {
        if seen
            .iter()
            .any(|m: &SymmetricBitMatrix<N>| m.is_equivalent_to(&mat, hash))
        {
            continue;
        }
        seen.push(mat);
        unique_matrices.push(mat);
    }
    unique_matrices
}
