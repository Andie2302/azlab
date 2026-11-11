//! Lineare Algebra fester Größe.
//!
//! Vektoren und Matrizen liegen auf dem Stack, ihre Größe steht im Typ.
//! Nichts in diesem Modul braucht Heap-Speicher.

mod error;
mod vector;

pub use error::LinalgError;
pub use vector::Vector;
