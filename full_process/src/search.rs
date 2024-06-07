use fxhash::{FxHashMap, FxHashSet};
use itertools::Itertools;

use crate::matrix::{MatrixHash, SymmetricBitMatrix};

#[derive(Debug, Clone)]
pub struct RowOrderStore<const N: usize> {
    memory: FxHashMap<[u64; N], Vec<[usize; N]>>,
}

impl<const N: usize> Default for RowOrderStore<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> RowOrderStore<N> {
    pub fn new() -> Self {
        let mut memory = FxHashMap::default();
        let mut hash_array = [0; N];
        for deltas in 0..(1 << (N - 1)) {
            for n in 0..(N - 1) {
                hash_array[n + 1] = hash_array[n] + (deltas >> n & 1);
            }
            memory.insert(hash_array, Self::generate(&hash_array));
        }
        Self { memory }
    }

    fn generate(hash: &[u64]) -> Vec<[usize; N]> {
        let mut row_orders = Vec::new();
        Self::generate_impl(hash, &mut [0; N], &mut row_orders);
        row_orders
    }

    fn generate_impl(hash: &[u64], row_order: &mut [usize; N], row_orders: &mut Vec<[usize; N]>) {
        if hash.is_empty() {
            row_orders.push(*row_order);
            return;
        }

        let same_count = hash.iter().take_while(|&&x| x == hash[0]).count();
        let residual = &hash[same_count..];
        let i0 = N - hash.len();
        for perm in (i0..(i0 + same_count)).permutations(same_count) {
            for (i, &p) in perm.iter().enumerate() {
                row_order[i0 + i] = p;
            }
            Self::generate_impl(residual, row_order, row_orders);
        }
    }

    pub fn get(&self, hash_array: &[u64; N]) -> Option<&Vec<[usize; N]>> {
        self.memory.get(hash_array)
    }
}

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

pub fn make_unique<const N: usize>(
    matrices: &[SymmetricBitMatrix<N>],
    hash: &MatrixHash<N>,
    store: &RowOrderStore<N>,
) -> Vec<SymmetricBitMatrix<N>> {
    let mut unique_matrices = Vec::new();
    let mut seen = Vec::new();
    for mat in matrices {
        if seen
            .iter()
            .any(|s: &FxHashSet<SymmetricBitMatrix<N>>| s.contains(mat))
        {
            continue;
        }

        let variants = make_variants(mat, hash, store);
        seen.push(variants);
        unique_matrices.push(*mat);
    }
    unique_matrices
}

fn make_variants<const N: usize>(
    matrix: &SymmetricBitMatrix<N>,
    hash: &MatrixHash<N>,
    store: &RowOrderStore<N>,
) -> FxHashSet<SymmetricBitMatrix<N>> {
    hash.generate_row_orders(store)
        .iter()
        .map(|row_order| matrix.create_rearranged(row_order))
        .collect()
}
