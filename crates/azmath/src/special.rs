//! Spezielle Funktionen: Gamma-, Beta- und Fehlerfunktion.
//!
//! | Funktion | Bedeutung |
//! |----------|-----------|
//! | [`gamma`], [`ln_gamma`] | Gammafunktion `Γ(x)` und ihr Logarithmus |
//! | [`beta`] | Betafunktion `B(a, b)` |
//! | [`gamma_p`], [`gamma_q`] | regularisierte unvollständige Gammafunktionen |
//! | [`erf`], [`erfc`] | Fehlerfunktion und ihr Komplement |

mod error_function;
mod gamma;

#[cfg(test)]
mod tests;

pub use error_function::{erf, erfc};
pub use gamma::{beta, gamma, gamma_p, gamma_q, ln_gamma};
