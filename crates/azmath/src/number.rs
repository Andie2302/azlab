//! Zahlentypen oberhalb der primitiven Gleitkommazahlen.
//!
//! Alle Typen implementieren die Traits aus [`crate::scalar`] und lassen sich
//! deshalb in jeden generischen Algorithmus einsetzen. Nichts in diesem Modul
//! braucht Heap-Speicher.

#[macro_use]
mod macros;

mod cayley_dickson;
mod complex;
mod dual;

pub use cayley_dickson::{CayleyDickson, Octonion};
pub use complex::Complex;
pub use dual::{Dual, DualN};
