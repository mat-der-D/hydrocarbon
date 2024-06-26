use std::hash::{Hash, Hasher};

use fxhash::FxHasher;

use crate::ordering::RowOrderStore;

#[derive(Debug, Clone, Copy)]
pub struct SymmetricBitMatrix<const N: usize> {
    rows: [u16; N],
}

impl<const N: usize> SymmetricBitMatrix<N> {
    pub fn zero() -> Self {
        Self { rows: [0; N] }
    }

    fn bit_at(&self, row: usize, col: usize) -> u16 {
        self.rows[row] >> col & 1
    }

    fn flip(&mut self, row: usize, col: usize) {
        self.rows[row] ^= 1 << col;
        self.rows[col] ^= 1 << row;
    }

    const REPUNIT: u16 = {
        let mut repunit = 0;
        let mut i = 0;
        while i < N {
            repunit = repunit << 1 | 1;
            i += 1;
        }
        repunit
    };

    fn is_connected(&self) -> bool {
        let mut connected = 1;
        for _ in 0..N {
            let tmp = connected;
            for (i, &row) in self.rows.iter().enumerate() {
                if tmp & 1 << i != 0 {
                    connected |= row;
                }
            }
            if connected == Self::REPUNIT {
                return true;
            }
        }
        false
    }

    pub fn create_rearranged(&self, row_order: &[usize; N]) -> Self {
        let mut rows_new = [0; N];
        for (row_new, &i_old) in rows_new.iter_mut().zip(row_order.iter()) {
            let row_old = unsafe { self.rows.get_unchecked(i_old) };
            for (j_new, &j_old) in row_order.iter().enumerate() {
                *row_new |= (row_old >> j_old & 1) << j_new;
            }
        }
        Self { rows: rows_new }
    }

    const UNIT_MATRIX: [[u16; N]; N] = {
        let mut mat = [[0; N]; N];
        let mut i = 0;
        while i < N {
            mat[i][i] = 1;
            i += 1;
        }
        mat
    };

    fn calc_features(&self) -> [u64; N] {
        const STEPS: usize = 3;
        let mut raw_features = [[0; STEPS]; N];
        let mut mat = Self::UNIT_MATRIX;
        for s in 0..N.min(STEPS) {
            let mut new_mat = [[0; N]; N];

            for (new_row, row_bits) in new_mat.iter_mut().zip(self.rows.iter()) {
                for (j, row) in mat.iter().enumerate() {
                    if row_bits & 1 << j == 0 {
                        continue;
                    }
                    for (new_elem, elem) in new_row.iter_mut().zip(row.iter()) {
                        *new_elem += elem;
                    }
                }
            }

            for (i, (raw_feat, new_row)) in raw_features.iter_mut().zip(new_mat.iter()).enumerate()
            {
                raw_feat[s] = {
                    let sqsum: u16 = new_row.iter().map(|&x| x * x).sum();
                    (sqsum as u32) ^ (new_row[i] as u32) << 16
                };
            }
            mat = new_mat;
        }

        let mut features = [0; N];
        for (feat, raw_feat) in features.iter_mut().zip(raw_features.iter()) {
            let mut hasher = FxHasher::default();
            raw_feat.hash(&mut hasher);
            *feat = hasher.finish();
        }
        features
    }

    const INDEX_ARRAY: [usize; N] = {
        let mut array = [0; N];
        let mut i = 0;
        while i < N {
            array[i] = i;
            i += 1;
        }
        array
    };

    pub fn canonicalize(&self) -> (Self, MatrixHash<N>) {
        let features = self.calc_features();
        let mut row_order = Self::INDEX_ARRAY;
        row_order.sort_by_key(|&i| features[i]);
        let canon = self.create_rearranged(&row_order);
        let mut sorted_features = [0; N];
        for (hash, &i_old) in sorted_features.iter_mut().zip(row_order.iter()) {
            *hash = features[i_old];
        }
        (canon, MatrixHash::new(sorted_features))
    }
}

impl<const N: usize> std::fmt::Display for SymmetricBitMatrix<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for row in self.rows.iter() {
            for i in 0..N {
                write!(f, "{}", row >> i & 1)?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct HydroCarbonMatrixIter<const N: usize> {
    current: SymmetricBitMatrix<N>,
    row_hashes: [u32; N],
    row: usize,
    col: usize,
    min_row: usize,
    min_col: usize,
}

impl<const N: usize> Default for HydroCarbonMatrixIter<N> {
    fn default() -> Self {
        Self::try_from_matrix(SymmetricBitMatrix::zero(), 0, 1).unwrap()
    }
}

impl<const N: usize> HydroCarbonMatrixIter<N> {
    pub fn create_multi(digits: usize) -> Vec<Self> {
        fn _next_row_col<const N: usize>(row: usize, col: usize) -> (usize, usize) {
            if col == N - 1 {
                (row + 1, row + 2)
            } else {
                (row, col + 1)
            }
        }

        let digits = digits.min(N * (N - 1) / 2 - 1);
        let n_max = 1 << digits;
        let mut iters = Vec::with_capacity(n_max);

        for n in 0..n_max {
            let mut mat = SymmetricBitMatrix::zero();
            let mut row = 0;
            let mut col = 1;
            for i in 0..digits {
                if n & 1 << i != 0 {
                    mat.flip(row, col);
                }
                (row, col) = _next_row_col::<N>(row, col);
            }
            if let Ok(iter) = Self::try_from_matrix(mat, row, col) {
                iters.push(iter);
            }
        }

        iters
    }

    fn try_from_matrix(
        mat: SymmetricBitMatrix<N>,
        min_row: usize,
        min_col: usize,
    ) -> Result<Self, &'static str> {
        let mut row_hashes = [0; N];
        for (hash, &row) in row_hashes.iter_mut().zip(mat.rows.iter()).take(min_row) {
            *hash = (row.count_ones() << 16) | (row as u32);
        }
        for (i, &row) in row_hashes.iter().enumerate().take(min_row) {
            if i == 0 && row > 0x4f000 || i != 0 && row > row_hashes[i - 1] {
                return Err("Invalid matrix");
            }
        }

        Ok(Self {
            current: mat,
            row_hashes,
            row: min_row,
            col: min_col,
            min_row,
            min_col,
        })
    }

    fn is_fine(&mut self) -> bool {
        if self.col != N - 1 && self.current.bit_at(self.row, self.col) == 0 {
            return true;
        }

        let next_row = self.row + (self.col == N - 1) as usize;
        let mut prev = if self.row == 0 {
            0x4f000
        } else {
            self.row_hashes[self.row - 1]
        };
        for (i, &row) in self.current.rows.iter().enumerate().skip(self.row) {
            let this = (row.count_ones() << 16) | (row as u32);
            if this > prev {
                return false;
            }
            if i < next_row {
                prev = this;
                self.row_hashes[i] = this;
            }
        }
        self.row < N - 2 || self.current.is_connected()
    }

    fn step_forward(&mut self) -> bool {
        if self.row >= N - 2 {
            return false;
        }
        if self.col == N - 1 {
            self.row += 1;
            self.col = self.row + 1;
        } else {
            self.col += 1;
        }
        true
    }

    fn step_backward(&mut self) -> bool {
        if (self.row, self.col) == (self.min_row, self.min_col) {
            return false;
        }
        if self.col == self.row + 1 {
            self.row -= 1;
            self.col = N - 1;
        } else {
            self.col -= 1;
        }
        true
    }
}

impl<const N: usize> PartialEq for SymmetricBitMatrix<N> {
    fn eq(&self, other: &Self) -> bool {
        if N <= 11 {
            u64::from(*self) == u64::from(*other)
        } else {
            u128::from(*self) == u128::from(*other)
        }
    }
}

impl<const N: usize> Eq for SymmetricBitMatrix<N> {}

macro_rules! impl_from_matrix_into_hash {
    ($num_type:ty) => {
        impl<const N: usize> From<SymmetricBitMatrix<N>> for $num_type {
            fn from(mat: SymmetricBitMatrix<N>) -> Self {
                let mut hash = 0;
                let mut shift = 0;
                let mut mask = 0;
                for (i, &row) in mat.rows.iter().enumerate() {
                    hash |= ((row & mask) as $num_type) << shift;
                    shift += i;
                    mask = mask << 1 | 1;
                }
                hash
            }
        }
    };
}

impl_from_matrix_into_hash!(u64); // N <= 11 の場合のみハッシュが重複しないことが保証される
impl_from_matrix_into_hash!(u128); // N <= 16 で常にハッシュが重複しないことが保証される

impl<const N: usize> Iterator for HydroCarbonMatrixIter<N> {
    type Item = SymmetricBitMatrix<N>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut skip_judgement = self.row == N - 2;
        let mut is_forward = true;
        loop {
            if !skip_judgement && is_forward && self.is_fine() {
                let result = self.step_forward();
                if !result {
                    break;
                }
            } else {
                self.current.flip(self.row, self.col);
                if self.current.bit_at(self.row, self.col) == 0 {
                    let result = self.step_backward();
                    if !result {
                        return None;
                    }
                    is_forward = false;
                } else {
                    is_forward = true;
                }
            }
            skip_judgement = false;
        }
        Some(self.current)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct MatrixHash<const N: usize> {
    features: [u64; N],
    canonical_features: [u64; N],
}

impl<const N: usize> MatrixHash<N> {
    fn new(features: [u64; N]) -> Self {
        let canonical_features = Self::canonicalize(&features);
        Self {
            features,
            canonical_features,
        }
    }

    fn canonicalize(hash_array: &[u64; N]) -> [u64; N] {
        let mut canonicalized = [0; N];
        let mut n = 0;
        let mut memory = hash_array[0];
        for (canon, &h) in canonicalized.iter_mut().zip(hash_array.iter()) {
            if h != memory {
                n += 1;
                memory = h;
            }
            *canon = n;
        }
        canonicalized
    }

    pub fn generate_row_orders<'a>(&'a self, store: &'a RowOrderStore<N>) -> &Vec<[usize; N]> {
        store.get(&self.canonical_features).unwrap()
    }
}

impl<const N: usize> PartialEq for MatrixHash<N> {
    fn eq(&self, other: &Self) -> bool {
        self.features == other.features
    }
}

impl<const N: usize> Eq for MatrixHash<N> {}

impl<const N: usize> Hash for MatrixHash<N> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.features.hash(state);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SymmetricTwoBitsMatrix<const N: usize> {
    rows: [u32; N],
}

impl<const N: usize> SymmetricTwoBitsMatrix<N> {
    pub fn create_rearranged(&self, row_order: &[usize]) -> Self {
        let mut rows_new = [0; N];
        for (row_new, i_old) in rows_new.iter_mut().zip(row_order.iter()) {
            let row_old = unsafe { self.rows.get_unchecked(*i_old) };
            for (j_new, j_old) in row_order.iter().enumerate() {
                *row_new |= (row_old >> (2 * j_old) & 0b11) << (2 * j_new);
            }
        }
        Self { rows: rows_new }
    }

    fn count_bonds(row: u32) -> u32 {
        2 * (row & 0xaaaa_aaaa).count_ones() + (row & 0x5555_5555).count_ones()
    }

    pub fn count_hydrogen(&self) -> u32 {
        let num_c: u32 = self.rows.iter().map(|&row| Self::count_bonds(row)).sum();
        4 * N as u32 - num_c
    }

    pub fn create_dehydrogenated_unchecked(&self, row: usize, col: usize) -> Self {
        let mut new_matrix = *self;
        new_matrix.rows[row] += 1 << (2 * col);
        new_matrix.rows[col] += 1 << (2 * row);
        new_matrix
    }

    const MAX_BONDS: u32 = if N == 2 { 3 } else { 4 };
    pub fn dehydrogenatable_bonds(&self) -> Vec<(usize, usize)> {
        let mut is_hydrogenatable = [false; N];
        let mut bonds = Vec::with_capacity(2 * N);
        for (i, &row) in self.rows.iter().enumerate() {
            if Self::count_bonds(row) == Self::MAX_BONDS {
                continue;
            }
            is_hydrogenatable[i] = true;
            for (j, able) in is_hydrogenatable.iter().enumerate().take(i) {
                if *able && row & (0b11 << (2 * j)) != 0 {
                    bonds.push((i, j));
                }
            }
        }
        bonds
    }
}

impl<const N: usize> From<SymmetricBitMatrix<N>> for SymmetricTwoBitsMatrix<N> {
    fn from(mat: SymmetricBitMatrix<N>) -> Self {
        let mut rows = [0; N];
        for (&row, row_new) in mat.rows.iter().zip(rows.iter_mut()) {
            let row_u32 = row as u32;
            for i in 0..N {
                *row_new |= (row_u32 << i) & (1 << (2 * i));
            }
        }
        Self { rows }
    }
}
