//! Lineare Algebra.
//!
//! - Feste Größe: [`Vector`] und [`Matrix`] liegen auf dem Stack, ihre Größe
//!   steht im Typ.
//! - Größe zur Laufzeit: [`MatrixView`], [`MatrixViewMut`] und [`CsrMatrix`]
//!   sind Sichten auf Puffer des Aufrufers.
//! - [`LinearOperator`] verbindet alle Speicherformen mit den iterativen
//!   Lösern.
//!
//! Nichts in diesem Modul braucht Heap-Speicher.

mod decomp;
mod error;
mod matrix;
mod operator;
mod sparse;
mod vector;
mod view;

#[cfg(test)]
mod dynamic_tests;

pub use decomp::{Cholesky, Lu, Qr};
pub use error::LinalgError;
pub use matrix::Matrix;
pub use operator::LinearOperator;
pub use sparse::CsrMatrix;
pub use vector::Vector;
pub use view::{MatrixView, MatrixViewMut};
