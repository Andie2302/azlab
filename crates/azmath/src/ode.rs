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
//! Für Systeme, deren Größe erst zur Laufzeit feststeht (etwa aus der
//! Ortsdiskretisierung einer partiellen Differentialgleichung), gibt es die
//! Einzelschritte [`euler_slice`] und [`rk4_slice`] auf Slices.

mod adaptive;
mod error;
mod fixed;
mod implicit;
mod slices;

#[cfg(test)]
mod tests;

pub use adaptive::{OdeOptions, dormand_prince, dormand_prince_observed};
pub use error::{OdeError, OdeResult, OdeState};
pub use fixed::{euler, euler_step, rk4, rk4_step};
pub use implicit::{implicit_euler, trapezoid};
pub use slices::{euler_slice, rk4_slice};
