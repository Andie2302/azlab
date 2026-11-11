//! Cayley-Dickson-Konstruktion: Verdopplung einer Algebra.
//!
//! Aus den reellen Zahlen entstehen so nacheinander die komplexen Zahlen,
//! die Quaternionen, die Oktonionen und so weiter. Mit jeder Stufe geht eine
//! Eigenschaft verloren, und die Trait-Implementierungen bilden das ab:
//!
//! | Stufe | Beispiel | Traits |
//! |-------|----------|--------|
//! | 1 | `CayleyDickson<f64>` (komplex) | `Field`, `Commutative` |
//! | 2 | `CayleyDickson<CayleyDickson<f64>>` (Quaternion) | `Field`, nicht kommutativ |
//! | 3 | [`Octonion`] | nur `Algebra`, nicht assoziativ |
//!
//! Für den täglichen Gebrauch gibt es die handlicheren Typen
//! [`Complex`](crate::number::Complex) und
//! [`Quaternion`](crate::number::Quaternion) mit Umwandlungen.

mod convert;
mod ops;
mod repr;
mod traits;

#[cfg(test)]
mod tests;

pub use repr::{CayleyDickson, Octonion};
