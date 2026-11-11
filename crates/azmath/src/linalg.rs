//! Lineare Algebra fester Größe.
//!
//! Vektoren und Matrizen liegen auf dem Stack, ihre Größe steht im Typ.
//! Nichts in diesem Modul braucht Heap-Speicher.

mod decomp;
mod error;
mod matrix;
mod vector;

pub use decomp::{Cholesky, Lu, Qr};
pub use error::LinalgError;
pub use matrix::Matrix;
pub use vector::Vector;
