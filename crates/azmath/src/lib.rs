//! # azmath
//!
//! Reine Mathematik-Bibliothek des azlab-Workspace: Zahlentypen,
//! Rechenmethoden und Löser. Fachliche Crates (Stoffdaten, Strömung,
//! Robotik, KI, ...) bauen darauf auf, azmath selbst kennt keine Fachdomäne.
//!
//! ## Grundsätze
//!
//! - `no_std` ohne `alloc` ist der Normalfall. Alles, was Heap oder `std`
//!   braucht, liegt hinter den Features `alloc` bzw. `std`.
//! - Algorithmen werden generisch über die Traits aus [`scalar`] geschrieben,
//!   nie direkt gegen `f32`/`f64`. Neue Zahlentypen (Dualzahlen, komplexe
//!   Zahlen, später `f16`/`f128`) funktionieren dadurch ohne Änderung an
//!   bestehendem Code.
//! - Wer die Elementarfunktionen tatsächlich berechnet (`libm`, `std`, ...),
//!   ist ein Implementierungsdetail des privaten Backend-Moduls und wird nur
//!   über Feature-Flags gewählt.
//!
//! ## Features
//!
//! | Feature | Wirkung |
//! |---------|---------|
//! | `libm` (Standard) | Elementarfunktionen über die Crate `libm`, läuft in `no_std` |
//! | `std`   | Elementarfunktionen über `std`, aktiviert auch `alloc` |
//! | `alloc` | Heap-basierte Typen (dynamische Größen) |

#![no_std]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

#[macro_use]
mod macros;

pub mod autodiff;
mod backend;
pub mod interp;
pub mod interval;
pub mod lie;
pub mod linalg;
pub mod number;
pub mod ode;
pub mod poly;
pub mod quad;
pub mod scalar;
pub mod solve;
pub mod tensor;

#[cfg(test)]
mod testutil;

pub use linalg::{
    Cholesky, CsrMatrix, LinalgError, LinearOperator, Lu, Matrix, MatrixView, MatrixViewMut, Qr,
    Vector,
};
pub use number::{CayleyDickson, Complex, Dual, DualN, Octonion, Quaternion};
pub use interval::Interval;
pub use poly::Polynomial;
pub use scalar::{Algebra, Commutative, Conjugate, Elementary, Field, Real, RealElementary, Ring};
pub use solve::{Solution, SolveError, Tolerance};
pub use tensor::{TensorView, TensorViewMut};
