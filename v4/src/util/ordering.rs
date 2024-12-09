use fxhash::FxHashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Permutation<const N: usize> {
    ordering: [usize; N],
}

impl<const N: usize> Permutation<N> {
    pub const fn new(ordering: [usize; N]) -> Self {
        Self { ordering }
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

    pub const IDENTITY: Self = Self::new(Self::INDEX_ARRAY);

    fn new_transposition(i: usize, j: usize) -> Self {
        let mut ordering = Self::INDEX_ARRAY;
        ordering.swap(i, j);
        Self::new(ordering)
    }

    pub fn ordering(&self) -> &[usize; N] {
        &self.ordering
    }

    pub fn inverse(&self) -> Self {
        let mut inverse = [0; N];
        for i in 0..N {
            inverse[self.ordering[i]] = i;
        }
        Self::new(inverse)
    }
}

impl<const N: usize> std::ops::Mul for Permutation<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut ordering = [0; N];
        for (item, &rhs_item) in ordering.iter_mut().zip(rhs.ordering.iter()) {
            *item = self.ordering[rhs_item];
        }
        Self::new(ordering)
    }
}

pub trait Permutable<const N: usize> {
    fn permute_by(&self, permutation: &Permutation<N>) -> Self;
}

#[derive(Debug, Clone)]
pub struct PermutationsStore<const N: usize> {
    memory: FxHashMap<[u64; N], Vec<Permutation<N>>>,
}

impl<const N: usize> Default for PermutationsStore<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> PermutationsStore<N> {
    pub fn new() -> Self {
        let mut memory = FxHashMap::default();
        for deltas in 0..(1 << (N - 1)) {
            let hash_array = Self::deltas_to_hash_array(deltas);
            memory.insert(hash_array, Self::generate(&hash_array));
        }
        Self { memory }
    }

    fn deltas_to_hash_array(deltas: u64) -> [u64; N] {
        let mut hash_array = [0; N];
        for n in 0..(N - 1) {
            hash_array[n + 1] = hash_array[n] + (deltas >> n & 1);
        }
        hash_array
    }

    fn generate(hash: &[u64; N]) -> Vec<Permutation<N>> {
        let mut perms = Vec::new();
        let mut i_start = 0;
        for i_end in 1..N {
            if hash[i_start] != hash[i_end] {
                i_start = i_end;
            } else {
                perms.push(Permutation::new_transposition(i_start, i_end));
            }
        }
        perms
    }

    pub fn get(&self, key: &[u64; N]) -> &[Permutation<N>] {
        &self.memory[key]
    }
}
