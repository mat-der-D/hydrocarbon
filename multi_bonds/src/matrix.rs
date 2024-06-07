use std::{fmt, fs, io, path::Path, str::FromStr};

use serde::Deserialize;

use crate::search::RowOrderStore;

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
}

impl<const N: usize> FromStr for SymmetricTwoBitsMatrix<N> {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() != N * N {
            return Err("invalid length".to_string());
        }

        let mut rows = [0; N];
        let mut chars = s.chars();
        for row in rows.iter_mut() {
            for j in 0..N {
                *row |= match chars.next().unwrap() {
                    '0' => 0b00,
                    '1' => 0b01,
                    '2' => 0b10,
                    '3' => 0b11,
                    _ => return Err("invalid character".to_string()),
                } << (2 * j);
            }
        }
        Ok(Self { rows })
    }
}

impl<const N: usize> fmt::Display for SymmetricTwoBitsMatrix<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for row in self.rows.iter() {
            for j in 0..N {
                let bits = (row >> (2 * j)) & 0b11;
                let c = match bits {
                    0b00 => '0',
                    0b01 => '1',
                    0b10 => '2',
                    0b11 => '3',
                    _ => unreachable!(),
                };
                write!(f, "{}", c)?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl<'de, const N: usize> Deserialize<'de> for SymmetricTwoBitsMatrix<N> {
    fn deserialize<D>(deserializer: D) -> Result<SymmetricTwoBitsMatrix<N>, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        Self::from_str(&s).map_err(serde::de::Error::custom)
    }
}

#[derive(Debug, Clone)]
struct MatrixAndHash<const N: usize> {
    matrix: SymmetricTwoBitsMatrix<N>,
    hash: [u64; N],
}

impl<'de, const N: usize> Deserialize<'de> for MatrixAndHash<N> {
    fn deserialize<D>(deserializer: D) -> Result<MatrixAndHash<N>, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let value = serde_json::Value::deserialize(deserializer)?;
        let matrix = SymmetricTwoBitsMatrix::deserialize(value["matrix"].clone()).unwrap();
        let hash_vec = value["hash"].as_array().unwrap();
        let mut hash = [0; N];
        for (i, val) in hash_vec.iter().enumerate() {
            hash[i] = val.as_u64().unwrap();
        }
        Ok(MatrixAndHash { matrix, hash })
    }
}

pub fn load_matrices<const N: usize>(
    path: impl AsRef<Path>,
) -> io::Result<Vec<(SymmetricTwoBitsMatrix<N>, [u64; N])>> {
    let file = fs::File::open(&path)?;
    let reader = io::BufReader::new(file);
    let data: Vec<MatrixAndHash<N>> = serde_json::from_reader(reader)?;
    Ok(data.into_iter().map(|x| (x.matrix, x.hash)).collect())
}

pub fn get_symmetry<const N: usize>(
    matrix: &SymmetricTwoBitsMatrix<N>,
    hash: &[u64; N],
    store: &RowOrderStore<N>,
) -> Vec<[usize; N]> {
    let mut symmetry = Vec::new();
    let row_orders = store.get(hash).unwrap();
    for row_order in row_orders {
        let rearranged = matrix.create_rearranged(row_order);
        if rearranged == *matrix {
            symmetry.push(*row_order);
        }
    }
    symmetry
}
