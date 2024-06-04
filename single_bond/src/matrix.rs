use std::{
    fmt,
    hash::{Hash, Hasher},
};

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

impl<const N: usize> fmt::Display for SymmetricBitMatrix<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn _bit_at<const N: usize>(mat: &SymmetricBitMatrix<N>, row: usize, col: usize) -> u16 {
            if mat.rows[row] & (1 << col) != 0 {
                1
            } else {
                0
            }
        }

        for row in 0..N {
            for col in 0..N {
                write!(f, "{}", _bit_at(self, row, col))?;
            }
            writeln!(f)?;
        }
        Ok(())
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

impl<const N: usize> From<MatrixHash<N>> for u64 {
    fn from(hash: MatrixHash<N>) -> u64 {
        let mut ret = 0u64;
        for &h in hash.morgan_hashes.iter() {
            ret = ret.wrapping_add(h);
        }
        ret.wrapping_add(hash.eig3_hash)
    }
}
