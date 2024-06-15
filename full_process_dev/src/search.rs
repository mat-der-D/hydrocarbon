use std::hash::Hash;

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
            for (row_order_elem, &p) in row_order[i0..].iter_mut().zip(perm.iter()) {
                *row_order_elem = p;
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
        let mut prev = 0x4ffff;
        for (i, &row) in matrix
            .rows()
            .iter()
            .enumerate()
            .skip(idx_row.saturating_sub(2))
        {
            let this = (row.count_ones() << 16) | (row as u32);
            if this > prev {
                return false;
            }
            if i < idx_row {
                prev = this;
            }
        }
        idx_row < N - 1 || matrix.is_connected()
    }
}

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
