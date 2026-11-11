//! Lineare Algebra.
//!
//! - Feste Größe: [`Vector`] und [`Matrix`] liegen auf dem Stack, ihre Größe
//!   steht im Typ.
//! - Größe zur Laufzeit: [`MatrixView`], [`MatrixViewMut`] und [`CsrMatrix`]
//!   sind Sichten auf Puffer des Aufrufers.
//! - [`LinearOperator`] verbindet alle Speicherformen mit den iterativen
//!   Lösern.
//!
//! Ohne Features braucht nichts in diesem Modul Heap-Speicher. Mit dem
//! Feature `alloc` kommt die besitzende `DMatrix` dazu.

mod decomp;
#[cfg(feature = "alloc")]
mod dmatrix;
mod error;
mod matrix;
mod operator;
mod sparse;
mod vector;
mod view;

#[cfg(test)]
mod dynamic_tests;

pub use decomp::{Cholesky, Lu, Qr, Svd, SymmetricEigen};
#[cfg(feature = "alloc")]
pub use dmatrix::DMatrix;
pub use error::LinalgError;
pub use matrix::Matrix;
pub use operator::{AffineOperator, LinearOperator};
pub use sparse::CsrMatrix;
pub(crate) use vector::scaled_norm;
pub use vector::Vector;
pub use view::{MatrixView, MatrixViewMut};
