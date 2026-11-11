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

#[cfg(feature = "std")]
extern crate std;

mod backend;
pub mod number;
pub mod scalar;

#[cfg(test)]
mod testutil;

pub use number::Complex;
pub use scalar::{Algebra, Commutative, Conjugate, Elementary, Field, Real, RealElementary, Ring};
