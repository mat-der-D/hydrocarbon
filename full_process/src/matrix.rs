use std::hash::{Hash, Hasher};

use fxhash::FxHasher;

use crate::search::RowOrderStore;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
            for i in 0..N {
                if tmp & 1 << i != 0 {
                    connected |= self.rows[i];
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
        for (row_new, i_old) in rows_new.iter_mut().zip(row_order.iter()) {
            let row_old = unsafe { self.rows.get_unchecked(*i_old) };
            for (j_new, j_old) in row_order.iter().enumerate() {
                *row_new |= (row_old >> j_old & 1) << j_new;
            }
        }
        Self { rows: rows_new }
    }

    pub fn make_hash(&self) -> MatrixHash<N> {
        MatrixHash::new(self.calc_morgan_hashes(), self.calc_eig3_hash())
    }

    fn calc_morgan_hashes(&self) -> [u64; N] {
        let mut hash_array = [[0; N]; N];
        let mut hash = [1u32; N];
        for s in 0..N {
            let mut new_hash = [0; N];
            for (i, new_hash_i) in new_hash.iter_mut().enumerate() {
                for (j, hash_j) in hash.iter().enumerate() {
                    if self.rows[i] & 1 << j != 0 {
                        *new_hash_i += hash_j;
                    }
                    hash_array[i][s] = *new_hash_i;
                }
            }
            hash = new_hash;
        }

        let mut ret = [0; N];
        for (ret_elem, arr) in ret.iter_mut().zip(hash_array.iter()) {
            let mut hasher = FxHasher::default();
            arr.hash(&mut hasher);
            *ret_elem = hasher.finish();
        }
        ret
    }

    fn calc_eig3_hash(&self) -> u64 {
        let mut eig3 = 0;
        for (i, row_i) in self.rows.iter().enumerate() {
            for j in 0..i {
                if row_i & (1 << j) == 0 {
                    continue;
                }
                for row_k in self.rows.iter().take(j) {
                    eig3 += (row_k >> i) & (row_k >> j) & 1;
                }
            }
        }

        let mut state = FxHasher::default();
        eig3.hash(&mut state);
        state.finish()
    }

    pub fn partially_canonicalize(&self) -> Self {
        let mut row_order = [0; N];
        for (i, row_order_i) in row_order.iter_mut().enumerate() {
            *row_order_i = i;
        }
        let morgan_hashes = self.calc_morgan_hashes();
        row_order.sort_by_key(|&i| morgan_hashes[i]);
        self.create_rearranged(&row_order)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MatrixHash<const N: usize> {
    morgan_hashes: [u64; N],
    canonical_morgan_hashes: [u64; N],
    eig3_hash: u64,
}

impl<const N: usize> MatrixHash<N> {
    fn new(morgan_hashes: [u64; N], eig3_hash: u64) -> Self {
        let canonical_morgan_hashes = Self::canonicalize_hashes(&morgan_hashes);
        Self {
            morgan_hashes,
            canonical_morgan_hashes,
            eig3_hash,
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
        let mut num_h = 0;
        for mut row in self.rows {
            for _ in 0..N {
                num_h += row & 0b11;
                row >>= 2;
            }
        }
        4 * N as u32 - num_h
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

    pub fn generate_symmetry(
        &self,
        hash: &MatrixHash<N>,
        store: &RowOrderStore<N>,
    ) -> Vec<[usize; N]> {
        let mut symmetry = Vec::new();
        let row_orders = hash.generate_row_orders(store);
        for row_order in row_orders {
            let rearranged = self.create_rearranged(row_order);
            if rearranged == *self {
                symmetry.push(*row_order);
            }
        }
        symmetry
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
