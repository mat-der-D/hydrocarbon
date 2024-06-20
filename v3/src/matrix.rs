use std::hash::{Hash, Hasher};

use fxhash::FxHasher;

use crate::ordering::RowOrderStore;

#[derive(Debug, Clone, Copy)]
pub struct SymmetricBitMatrix<const N: usize> {
    rows: [u16; N],
}

impl<const N: usize> SymmetricBitMatrix<N> {
    fn zero() -> Self {
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
        let mut visited = 1;
        for _ in 0..N {
            let tmp = visited;
            for (i, &row) in self.rows.iter().enumerate() {
                visited |= row * (tmp >> i & 1);
            }
            if visited == Self::REPUNIT {
                return true;
            }
        }
        false
    }

    pub fn create_rearranged(&self, row_order: &[usize; N]) -> Self {
        let mut rows = [0; N];
        for (new_row, &i_old) in rows.iter_mut().zip(row_order.iter()) {
            let old_row = unsafe { self.rows.get_unchecked(i_old) };
            for (j_new, &j_old) in row_order.iter().enumerate() {
                *new_row |= (old_row >> j_old & 1) << j_new;
            }
        }
        Self { rows }
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
            let new_mat = self * &mat;

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

    pub fn canonicalize(&self) -> (Self, MatrixFeatures<N>) {
        let features = self.calc_features();
        let mut row_order = Self::INDEX_ARRAY;
        row_order.sort_by_key(|&i| features[i]);
        let canon = self.create_rearranged(&row_order);
        let mut sorted_features = [0; N];
        for (feat, &i_old) in sorted_features.iter_mut().zip(row_order.iter()) {
            *feat = features[i_old];
        }
        (canon, MatrixFeatures::new(sorted_features))
    }
}

impl<const N: usize> std::ops::Mul<&[[u16; N]; N]> for &SymmetricBitMatrix<N> {
    type Output = [[u16; N]; N];

    fn mul(self, rhs: &[[u16; N]; N]) -> Self::Output {
        let mut out = [[0; N]; N];
        for (out_row, &row_bits) in out.iter_mut().zip(self.rows.iter()) {
            for (j, row) in rhs.iter().enumerate() {
                if row_bits & 1 << j == 0 {
                    continue;
                }
                for (out_elem, &elem) in out_row.iter_mut().zip(row.iter()) {
                    *out_elem += elem;
                }
            }
        }
        out
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

#[derive(Debug, Clone, Copy)]
pub struct MatrixFeatures<const N: usize> {
    features: [u64; N],
    order_store_key: [u64; N],
}

impl<const N: usize> MatrixFeatures<N> {
    fn new(features: [u64; N]) -> Self {
        Self {
            features,
            order_store_key: Self::make_key(&features),
        }
    }

    fn make_key(features: &[u64; N]) -> [u64; N] {
        let mut key = [0; N];
        let mut n = 0;
        let mut memory = features[0];
        for (key_elem, &feat) in key.iter_mut().zip(features.iter()) {
            n += (feat != memory) as u64;
            memory = feat;
            *key_elem = n;
        }
        key
    }

    pub fn generate_row_orders<'a>(&'a self, store: &'a RowOrderStore<N>) -> &Vec<[usize; N]> {
        store.get(&self.order_store_key)
    }
}

impl<const N: usize> PartialEq for MatrixFeatures<N> {
    fn eq(&self, other: &Self) -> bool {
        self.features == other.features
    }
}

impl<const N: usize> Eq for MatrixFeatures<N> {}

impl<const N: usize> Hash for MatrixFeatures<N> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.features.hash(state);
    }
}

#[derive(Debug)]
pub struct HydroCarbonMatrixIter<const N: usize> {
    matrix: SymmetricBitMatrix<N>,
    row_hashes: [u32; N],
    row: usize,
    col: usize,
    min_row: usize,
    min_col: usize,
}

impl<const N: usize> HydroCarbonMatrixIter<N> {
    pub fn new_split(digits: usize) -> anyhow::Result<Vec<Self>> {
        if N == 0 {
            return Err(anyhow::anyhow!("N must be greater than 0"));
        }
        if N == 1 {
            todo!("yet to implement");
        }
        if digits > 32.min(N * (N - 1) / 2) {
            return Err(anyhow::anyhow!("Too many digits"));
        }

        let mut iters = Vec::new();

        let n_max = if digits == 32 {
            u32::MAX
        } else {
            1 << digits - 1
        };
        for n in 0..=n_max {
            let mut mat = SymmetricBitMatrix::zero();
            let mut row = 0;
            let mut col = 1;
            for i in 0..digits {
                if n & 1 << i != 0 {
                    mat.flip(row, col);
                }
                (row, col) = Self::next_row_col(row, col);
            }
            if let Ok(iter) = Self::try_from_matrix(mat, row, col) {
                iters.push(iter);
            }
        }

        Ok(iters)
    }

    fn try_from_matrix(
        matrix: SymmetricBitMatrix<N>,
        min_row: usize,
        min_col: usize,
    ) -> anyhow::Result<Self> {
        let mut row_hashes = [0; N];
        for (hash, &row) in row_hashes.iter_mut().zip(matrix.rows.iter()).take(min_row) {
            *hash = (row.count_ones() << 16) | (row as u32);
        }
        for (i, &row) in row_hashes.iter().enumerate().take(min_row) {
            if i == 0 && row > 0x4f000 || i != 0 && row > row_hashes[i - 1] {
                return Err(anyhow::anyhow!("Invalid row"));
            }
        }
        Ok(Self {
            matrix,
            row_hashes,
            row: min_row,
            col: min_col,
            min_row,
            min_col,
        })
    }

    fn next_row_col(row: usize, col: usize) -> (usize, usize) {
        if col < N - 1 {
            (row, col + 1)
        } else {
            (row + 1, row + 2)
        }
    }

    fn prev_row_col(row: usize, col: usize) -> (usize, usize) {
        match (row, col) {
            (0, 1) => (0, 0),
            (r, c) if r + 1 == c => (r - 1, N - 1),
            (r, c) => (r, c - 1),
        }
    }

    fn step_forward(&mut self) -> bool {
        (self.row, self.col) = Self::next_row_col(self.row, self.col);
        self.row < N - 1
    }

    fn step_backward(&mut self) -> bool {
        if (self.row, self.col) == (self.min_row, self.min_col) {
            false
        } else {
            (self.row, self.col) = Self::prev_row_col(self.row, self.col);
            true
        }
    }

    fn is_fine_up_to(&mut self, row: usize, col: usize) -> bool {
        if col < N - 1 && self.matrix.bit_at(row, col) == 0 {
            return true;
        }

        let (next_row, _) = Self::next_row_col(row, col);
        let mut prev = if row == 0 {
            0x4f000
        } else {
            self.row_hashes[row - 1]
        };

        for (i, &row) in self.matrix.rows.iter().enumerate().skip(row) {
            let this = (row.count_ones() << 16) | (row as u32);
            if this > prev {
                return false;
            }
            if i < next_row {
                prev = this;
                self.row_hashes[i] = this;
            }
        }
        next_row < N - 1 || self.matrix.is_connected()
    }

    const FINISHED_ROW_COL: (usize, usize) = (usize::MAX, usize::MAX);

    fn mark_as_finished(&mut self) {
        (self.row, self.col) = Self::FINISHED_ROW_COL;
    }

    fn is_finished(&self) -> bool {
        (self.row, self.col) == Self::FINISHED_ROW_COL
    }
}

impl<const N: usize> Iterator for HydroCarbonMatrixIter<N> {
    type Item = SymmetricBitMatrix<N>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_finished() {
            return None;
        }

        let mut is_forward = true;
        loop {
            if self.row >= N - 1 {
                if self.step_backward() {
                    is_forward = false;
                    continue;
                }

                let (prev_row, prev_col) = Self::prev_row_col(self.row, self.col);
                self.mark_as_finished();
                return if self.is_fine_up_to(prev_row, prev_col) {
                    Some(self.matrix)
                } else {
                    None
                };
            }

            if is_forward && self.is_fine_up_to(self.row, self.col) {
                if !self.step_forward() {
                    return Some(self.matrix);
                }
            } else {
                self.matrix.flip(self.row, self.col);
                if self.matrix.bit_at(self.row, self.col) == 0 {
                    if !self.step_backward() {
                        self.mark_as_finished();
                        return None;
                    }
                    is_forward = false;
                } else {
                    is_forward = true;
                }
            }
        }
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
    fn from(matrix: SymmetricBitMatrix<N>) -> Self {
        let mut rows = [0; N];
        for (new_row, &old_row) in rows.iter_mut().zip(matrix.rows.iter()) {
            let old_row_u32 = old_row as u32;
            for i in 0..N {
                *new_row |= (old_row_u32 << i) & (1 << (2 * i));
            }
        }
        Self { rows }
    }
}
