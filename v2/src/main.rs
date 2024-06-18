use fxhash::FxHashMap;
use matrix::HydroCarbonMatrixIter;

mod matrix;

fn main() {
    const N: usize = 10;
    let instant = std::time::Instant::now();
    let digits = 7;

    let mut handlers = Vec::new();
    let iters = HydroCarbonMatrixIter::<N>::create_multi(digits);
    for iter in iters {
        handlers.push(std::thread::spawn(move || {
            let mut hash2mats = FxHashMap::default();
            for mat in iter {
                let (mat, hash) = mat.canonicalize();
                hash2mats.entry(hash).or_insert_with(Vec::new).push(mat);
            }
            hash2mats
        }));
    }

    let mut hash2mats = FxHashMap::default();
    for handler in handlers {
        let h2m = handler.join().unwrap();
        for (hash, mut mats) in h2m {
            hash2mats
                .entry(hash)
                .or_insert_with(Vec::new)
                .append(&mut mats);
        }
    }

    println!("# of keys: {}", hash2mats.len());

    println!("{:?}", instant.elapsed());
}
