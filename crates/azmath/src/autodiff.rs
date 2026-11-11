//! Ableitungen ganzer Funktionen über Dualzahlen (Vorwärtsmodus).
//!
//! Die Funktion wird jeweils mit Dualzahlen statt gewöhnlichen Zahlen
//! aufgerufen. Man schreibt sie einmal generisch über die Skalar-Traits und
//! bekommt Wert und Ableitungen exakt bis auf Rundung, ohne Differenzenquotienten.

mod scalar;
mod vector;

#[cfg(test)]
mod tests;

pub use scalar::{derivative, second_derivative};
pub use vector::{gradient, hessian, jacobian, variables};
