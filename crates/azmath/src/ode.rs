//! Löser für Anfangswertprobleme gewöhnlicher Differentialgleichungen
//! `y' = f(t, y)`.
//!
//! | Verfahren | Art | Wofür |
//! |-----------|-----|-------|
//! | [`euler`], [`rk4`] | explizit, feste Schrittweite | einfache, nicht steife Probleme |
//! | [`dormand_prince`] | explizit, gesteuerte Schrittweite | Standard für nicht steife Probleme |
//! | [`implicit_euler`], [`trapezoid`] | implizit, feste Schrittweite | steife Probleme |
//!
//! Die Zahl der Gleichungen steht im Typ; nichts braucht Heap-Speicher.

mod adaptive;
mod error;
mod fixed;
mod implicit;

#[cfg(test)]
mod tests;

pub use adaptive::{OdeOptions, dormand_prince, dormand_prince_observed};
pub use error::{OdeError, OdeResult, OdeState};
pub use fixed::{euler, euler_step, rk4, rk4_step};
pub use implicit::{implicit_euler, trapezoid};
