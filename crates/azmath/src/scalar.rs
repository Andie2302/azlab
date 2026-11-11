//! Skalar-Traits: die algebraische Grundlage aller Algorithmen in azmath.
//!
//! Die Hierarchie ist absichtlich in kleine Stufen geteilt, damit jeder
//! Zahlentyp genau so viel zusagen muss, wie er leisten kann:
//!
//! | Trait | Zusage | Beispiele |
//! |-------|--------|-----------|
//! | [`Ring`] | `+ - *`, Null, Eins (Multiplikation muss nicht kommutativ sein) | Ganzzahlen, Quaternionen |
//! | [`Field`] | zusätzlich Division | `f32`, `f64`, komplexe Zahlen, Dualzahlen |
//! | [`Elementary`] | zusätzlich `sqrt`, `exp`, `sin`, ... (ohne Ordnung) | `f32`, `f64`, komplexe Zahlen, Dualzahlen |
//! | [`Real`] | Ordnung, Betrag, Konstanten | `f32`, `f64` |
//! | [`RealElementary`] | `atan2`, `hypot`, Runden, ... | `f32`, `f64` |
//!
//! Algorithmen fordern immer nur die schwächste Stufe, die sie brauchen.
//! `Ring`, `Field` und `Real` kommen ohne Rechen-Backend aus.

mod elementary;
mod field;
mod primitives;
mod real;
mod ring;

#[cfg(test)]
mod tests;

pub use elementary::{Elementary, RealElementary};
pub use field::Field;
pub use real::Real;
pub use ring::Ring;
