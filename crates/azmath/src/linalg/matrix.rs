//! Matrizen fester Größe.

mod expm;
mod ops;
mod repr;
mod square;

#[cfg(test)]
mod tests;

pub use repr::Matrix;
