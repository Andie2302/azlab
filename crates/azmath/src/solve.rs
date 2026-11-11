//! Löser für Gleichungen.
//!
//! - Eine Unbekannte: [`bisection()`], [`brent()`], [`secant()`], [`newton()`],
//!   [`newton_auto`], [`newton_bracketed`]. Keiner braucht ein Rechen-Backend.
//! - Mehrere Unbekannte: [`newton_system`], [`newton_system_with_jacobian`].
//!
//! Alle Löser nehmen ihre Abbruchkriterien als [`Tolerance`] und liefern eine
//! [`Solution`]. Konvergiert ein Verfahren nicht, ist das ein
//! [`SolveError`], der den letzten Stand enthält.

mod bisection;
mod brent;
mod error;
mod newton;
mod secant;
mod solution;
mod system;
mod tolerance;
mod util;

#[cfg(test)]
mod tests;

pub use bisection::bisection;
pub use brent::brent;
pub use error::{ScalarResult, SolveError};
pub use newton::{newton, newton_auto, newton_bracketed};
pub use secant::secant;
pub use solution::Solution;
pub use system::{SystemResult, newton_system, newton_system_with_jacobian};
pub use tolerance::Tolerance;
