//! Dualzahlen für automatisches Differenzieren im Vorwärtsmodus.
//!
//! Eine Dualzahl trägt neben dem Wert die Ableitung mit. Weil [`Dual`] und
//! [`DualN`] alle Skalar-Traits implementieren, liefert jeder generische
//! Algorithmus automatisch Ableitungen, wenn man ihn mit Dualzahlen aufruft.
//!
//! - [`Dual`]: eine Ableitungsrichtung.
//! - [`DualN`]: `N` Richtungen zugleich (Gradient), ohne Heap.
//! - `Dual<Dual<T>>`: zweite Ableitungen.
//!
//! Vergleiche (`==`, `<`, ...) berücksichtigen **nur den Wert**. Ein
//! Algorithmus verzweigt mit Dualzahlen deshalb genauso wie mit `f64`.

mod elementary;
mod ops;
mod repr;
mod tangent;
mod traits;

#[cfg(test)]
mod tests;

pub use repr::{Dual, DualN};
