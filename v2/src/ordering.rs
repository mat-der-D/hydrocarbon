use fxhash::FxHashMap;

#[derive(Debug, Clone)]
struct PermutationGenerator {
    current: Vec<usize>,
    is_first: bool,
}

impl PermutationGenerator {
    fn from_start_count(start: usize, count: usize) -> Self {
        let current = (start..(start + count)).collect();
        Self {
            current,
            is_first: true,
        }
    }
}

impl Iterator for PermutationGenerator {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.is_first {
            self.is_first = false;
            return Some(self.current.clone());
        }

        let mut i = self.current.len() - 1;
        while i > 0 && self.current[i - 1] >= self.current[i] {
            i -= 1;
        }
        if i == 0 {
            return None;
        }

        let mut j = self.current.len() - 1;
        while self.current[j] <= self.current[i - 1] {
            j -= 1;
        }

        self.current.swap(i - 1, j);
        self.current[i..].reverse();
        Some(self.current.clone())
    }
}

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
        for perm in PermutationGenerator::from_start_count(i0, same_count) {
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
