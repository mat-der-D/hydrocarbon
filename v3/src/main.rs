mod util;

use util::{display_counts, generate_hydrocarbons};

macro_rules! run_many {
    ($($n:expr),+;max_digits=$max_digits:expr;num_threads=$num_threads:expr) => {
        $(
            let digits = $max_digits.min($n * ($n - 1) / 2);
            let mat = generate_hydrocarbons::<$n>(digits, $num_threads)?;
            display_counts(&mat);
        )*
    };
}

fn main() -> anyhow::Result<()> {
    run_many!(1, 2, 3, 4, 5, 6, 7, 8; max_digits=5; num_threads=16);
    run_many!(9, 10; max_digits=7; num_threads=172);
    // run_many!(11; max_digits=7; num_threads=172);
    // run_many!(12; max_digits=7; num_threads=172);
    Ok(())
}
