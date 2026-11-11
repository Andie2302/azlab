//! Quaternionen mit Rotations-Grundfunktionen.
//!
//! Rotationsmatrizen und Euler-Winkel folgen mit den Matrix-Typen; ein
//! eigener Typ für Einheitsquaternionen kommt mit den Lie-Gruppen.

mod convert;
mod ops;
mod repr;
mod rotation;
mod traits;

#[cfg(test)]
mod tests;

pub use repr::Quaternion;
