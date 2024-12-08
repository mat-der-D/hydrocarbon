use fxhash::FxHashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Permutation<const N: usize> {
    ordering: [usize; N],
}

impl<const N: usize> Permutation<N> {
    pub fn new(ordering: [usize; N]) -> Self {
        Self { ordering }
    }

    pub fn identity() -> Self {
        let mut ordering = [0; N];
        for i in 0..N {
            ordering[i] = i;
        }
        Self::new(ordering)
    }

    pub fn is_identity(&self) -> bool {
        self.ordering.iter().enumerate().all(|(i, &x)| i == x)
    }

    pub fn new_transposition(i: usize, j: usize) -> Self {
        let mut ordering = [0; N];
        for k in 0..N {
            ordering[k] = k;
        }
        ordering.swap(i, j);
        Self::new(ordering)
    }

    pub fn ordering(&self) -> &[usize; N] {
        &self.ordering
    }

    pub fn permute<P>(&self, permutable: &P) -> P
    where
        P: Permutable<N>,
    {
        permutable.permute_by(self)
    }
}

impl<const N: usize> std::ops::Mul for Permutation<N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut ordering = [0; N];
        for i in 0..N {
            ordering[i] = self.ordering[rhs.ordering[i]];
        }
        Self::new(ordering)
    }
}

pub trait Permutable<const N: usize> {
    fn permute_by(&self, permutation: &Permutation<N>) -> Self;
}

#[derive(Debug, Clone)]
pub struct RowOrderStore<const N: usize> {
    memory: FxHashMap<[u64; N], Vec<Permutation<N>>>,
}

impl<const N: usize> Default for RowOrderStore<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> RowOrderStore<N> {
    pub fn new() -> Self {
        let mut memory = FxHashMap::default();
        for deltas in 0..(1 << (N - 1)) {
            let hash_array = Self::deltas_to_hash_array(deltas);
            memory.insert(hash_array, Self::generate(&hash_array));
        }
        Self { memory }
    }

    pub fn new_parallel(num_threads: u64) -> Self {
        let mut handlers = Vec::new();
        for i in 0..num_threads {
            let handler = std::thread::spawn(move || {
                let mut memory = FxHashMap::default();
                for deltas in (0..(1 << (N - 1))).filter(|x| x % num_threads == i) {
                    let hash_array = Self::deltas_to_hash_array(deltas);
                    memory.insert(hash_array, Self::generate(&hash_array));
                }
                memory
            });
            handlers.push(handler);
        }

        let mut memory = FxHashMap::default();

        for handler in handlers {
            let partial_memory = handler.join().unwrap();
            for (key, value) in partial_memory {
                memory.insert(key, value);
            }
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

    pub fn get(&self, key: &[u64; N]) -> &Vec<Permutation<N>> {
        &self.memory[key]
    }
}
