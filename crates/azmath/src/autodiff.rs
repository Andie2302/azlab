//! Automatisches Differenzieren.
//!
//! Die Funktion wird mit besonderen Zahlen statt gewöhnlichen aufgerufen.
//! Man schreibt sie einmal generisch über die Skalar-Traits und bekommt
//! Wert und Ableitungen exakt bis auf Rundung, ohne Differenzenquotienten.
//!
//! - **Vorwärtsmodus** über Dualzahlen: [`derivative`], [`gradient`],
//!   [`jacobian`], [`hessian`]. Günstig bei wenigen Eingängen.
//! - **Rückwärtsmodus** über ein Band ([`Tape`], [`Var`],
//!   [`reverse_gradient`]): ein Durchlauf liefert die Ableitungen nach
//!   allen Eingängen. Günstig bei vielen Eingängen und wenigen Ausgängen.

mod reverse;
mod scalar;
mod tape;
mod var;
mod vector;

#[cfg(test)]
mod reverse_tests;
#[cfg(test)]
mod tests;

pub use reverse::reverse_gradient;
pub use scalar::{derivative, second_derivative};
#[cfg(feature = "alloc")]
pub use tape::GrowTape;
pub use tape::{Record, Tape, TapeError};
pub use var::Var;
pub use vector::{gradient, hessian, jacobian, variables};
