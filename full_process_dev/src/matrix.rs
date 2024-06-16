use std::hash::{Hash, Hasher};

use fxhash::FxHasher;

use crate::search::RowOrderStore;

#[derive(Debug, Clone, Copy)]
pub struct SymmetricBitMatrix<const N: usize> {
    rows: [u16; N],
}

impl<const N: usize> SymmetricBitMatrix<N> {
    pub fn zero() -> Self {
        Self { rows: [0; N] }
    }

    pub fn rows(&self) -> &[u16; N] {
        &self.rows
    }

    const REPUNIT: u16 = {
        let mut repunit = 0;
        let mut i = 0;
        while i < N {
            repunit |= 1 << i;
            i += 1;
        }
        repunit
    };

    pub fn is_connected(&self) -> bool {
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

    pub fn flip(&mut self, row: usize, col: usize) {
        self.rows[row] ^= 1 << col;
        self.rows[col] ^= 1 << row;
    }

    pub fn is_symmetrical_under(&self, row_order: &[usize; N]) -> bool {
        for (&row_i_new, &i_old) in self.rows.iter().zip(row_order.iter()) {
            let &row_i_old = unsafe { self.rows.get_unchecked(i_old) };
            for (j_new, &j_old) in row_order.iter().enumerate() {
                if row_i_new >> j_new & 1 != row_i_old >> j_old & 1 {
                    return false;
                }
            }
        }
        true
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

        let mut ext_features = [0; N];
        for (ext_feat, &row) in ext_features.iter_mut().zip(self.rows.iter()) {
            *ext_feat = features.iter().enumerate().fold(0u64, |acc, (i, &f)| {
                acc.wrapping_add(f * ((row >> i) & 1) as u64)
            });
        }
        ext_features
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

    pub fn partially_canonicalize(&self) -> (Self, MatrixHash<N>) {
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
                for (i, &row) in mat.rows().iter().enumerate() {
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

    pub fn count_hydrogen(&self) -> u32 {
        let mut num_c = 0;
        for mut row in self.rows {
            for _ in 0..N {
                num_c += row & 0b11;
                row >>= 2;
            }
        }
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
        fn _count_bonds<const N: usize>(row: u32) -> u32 {
            2 * (row & 0xaaaa_aaaa).count_ones() + (row & 0x5555_5555).count_ones()
        }

        let mut is_hydrogenatable = [false; N];
        let mut bonds = Vec::with_capacity(2 * N);
        for (i, &row) in self.rows.iter().enumerate() {
            if _count_bonds::<N>(row) == Self::MAX_BONDS {
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
