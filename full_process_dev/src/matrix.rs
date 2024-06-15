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

    pub fn activate(&mut self, row: usize, col: usize) {
        self.rows[row] |= 1 << col;
        self.rows[col] |= 1 << row;
    }

    pub fn deactivate(&mut self, row: usize, col: usize) {
        self.rows[row] &= !(1 << col);
        self.rows[col] &= !(1 << row);
    }

    pub fn create_rearranged(&self, row_order: &[usize]) -> Self {
        let mut rows_new = [0; N];
        for (row_new, &i_old) in rows_new.iter_mut().zip(row_order.iter()) {
            let row_old = unsafe { self.rows.get_unchecked(i_old) };
            for (j_new, &j_old) in row_order.iter().enumerate() {
                *row_new |= (row_old >> j_old & 1) << j_new;
            }
        }
        Self { rows: rows_new }
    }

    pub fn make_hash(&self) -> MatrixHash<N> {
        let (morgan_hashes, traces_hash) = self.calc_hashes();
        MatrixHash::new(morgan_hashes, traces_hash)
    }

    fn calc_hashes(&self) -> ([u64; N], u64) {
        let mut morgans = [[0; N]; N];
        let mut traces = [0; N];

        // A^s を計算しながら, tr(A^s) と A^s [1, 1, ..., 1]^T を計算
        let mut mat = [[0; N]; N];
        for (i, row) in mat.iter_mut().enumerate() {
            row[i] = 1;
        }
        for (s, trace) in traces.iter_mut().enumerate() {
            let mut new_mat = [[0; N]; N];
            for (new_row, row) in new_mat.iter_mut().zip(mat.iter()) {
                for (new_row_elem, row_bits) in new_row.iter_mut().zip(self.rows.iter()) {
                    *new_row_elem = row
                        .iter()
                        .enumerate()
                        .map(|(j, elem)| elem * ((row_bits >> j) & 1))
                        .sum();
                }
            }
            for (morgan, row) in morgans.iter_mut().zip(new_mat.iter()) {
                morgan[s] = row.iter().map(|&elem| elem).sum();
            }
            *trace = new_mat.iter().enumerate().map(|(i, row)| row[i]).sum();
            mat = new_mat;
        }

        let mut hasher = FxHasher::default();
        traces.hash(&mut hasher);
        let traces_hash = hasher.finish();

        let mut morgan_hashes = [0; N];
        for (hash, morgan) in morgan_hashes.iter_mut().zip(morgans.iter()) {
            let mut hasher = FxHasher::default();
            morgan.hash(&mut hasher);
            *hash = hasher.finish();
        }
        (morgan_hashes, traces_hash)
    }

    pub fn partially_canonicalize(&self) -> (Self, MatrixHash<N>) {
        let mut row_order = [0; N];
        for (i, row_order_i) in row_order.iter_mut().enumerate() {
            *row_order_i = i;
        }
        let hash = self.make_hash();
        row_order.sort_by_key(|&i| hash.morgan_hashes[i]);
        let canon = self.create_rearranged(&row_order);
        let mut new_morgan_hahes = [0; N];
        for (new_hash, &i_old) in new_morgan_hahes.iter_mut().zip(row_order.iter()) {
            *new_hash = hash.morgan_hashes[i_old];
        }
        (canon, MatrixHash::new(new_morgan_hahes, hash.traces_hash))
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MatrixHash<const N: usize> {
    morgan_hashes: [u64; N],
    canonical_morgan_hashes: [u64; N],
    traces_hash: u64,
}

impl<const N: usize> MatrixHash<N> {
    fn new(morgan_hashes: [u64; N], traces_hash: u64) -> Self {
        let canonical_morgan_hashes = Self::canonicalize_hashes(&morgan_hashes);
        Self {
            morgan_hashes,
            canonical_morgan_hashes,
            traces_hash,
        }
    }

    fn canonicalize_hashes(hash_array: &[u64; N]) -> [u64; N] {
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
        store.get(&self.canonical_morgan_hashes).unwrap()
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
        fn _count_bonds<const N: usize>(mut row: u32) -> u32 {
            let mut num_bonds = 0;
            for _ in 0..N {
                num_bonds += row & 0b11;
                row >>= 2;
            }
            num_bonds
        }

        let mut is_hydrogenatable = [false; N];
        let mut bonds = Vec::new();
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
