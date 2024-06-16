use fxhash::FxHashMap;
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

    fn generate(hash: &[u64; N]) -> Vec<[usize; N]> {
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
        self.search_impl(
            0, 1, &mut f, 0x4f000, /* 0 行目がとりうる最大値 */
        );
    }

    fn search_impl(
        &mut self,
        row: usize,
        col: usize,
        f: &mut impl FnMut(&SymmetricBitMatrix<N>),
        conf: u32, // 最終確定行の情報を保持する. 前半はノードの次数, 後半は隣接ノードの情報
    ) {
        if row == N - 1 {
            f(&self.matrix);
            return;
        }

        let (next_row, next_col) = if col == N - 1 {
            (row + 1, row + 2)
        } else {
            (row, col + 1)
        };

        // 0 を入れる場合
        let mut conf_mut = conf;
        if next_row == row || Self::is_fine_up_to(&self.matrix, row, next_row, &mut conf_mut) {
            self.search_impl(next_row, next_col, f, conf_mut); // 次の最終確定行は直前の is_fine_up_to の結果を使う
        }
        self.matrix.flip(row, col);

        // 1 を入れる場合
        let mut conf_mut = conf;
        if Self::is_fine_up_to(&self.matrix, row, next_row, &mut conf_mut) {
            self.search_impl(next_row, next_col, f, conf_mut); // 次の最終確定行は直前の is_fine_up_to の結果を使う
        }
        self.matrix.flip(row, col);
    }

    fn is_fine_up_to(
        matrix: &SymmetricBitMatrix<N>,
        row: usize,
        next_row: usize,
        conf: &mut u32, // 最終確定行の情報を保持する. 前半はノードの次数, 後半は隣接ノードの情報
    ) -> bool {
        for (i, &row_bits) in matrix.rows().iter().enumerate().skip(row) {
            let this = (row_bits.count_ones() << 16) | (row_bits as u32);
            if this > *conf {
                return false;
            }
            if i < next_row {
                *conf = this;
            }
        }
        next_row < N - 1 || matrix.is_connected()
    }
}

pub fn make_symmetry<const N: usize>(
    mat: &SymmetricBitMatrix<N>,
    hash: &MatrixHash<N>,
    store: &RowOrderStore<N>,
) -> Vec<[usize; N]> {
    hash.generate_row_orders(store)
        .iter()
        .filter(|&&row_order| mat.is_symmetrical_under(&row_order))
        .copied()
        .collect()
}
