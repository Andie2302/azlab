//! Vektoren fester Länge.

mod ops;
mod products;
mod repr;

#[cfg(test)]
mod tests;

pub(crate) use products::scaled_norm;
pub use repr::Vector;
