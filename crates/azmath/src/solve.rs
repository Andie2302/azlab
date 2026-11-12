//! Löser für Gleichungen.
//!
//! - Eine Unbekannte: [`bisection()`], [`brent()`], [`secant()`], [`newton()`],
//!   [`newton_auto`], [`newton_bracketed`]. Keiner braucht ein Rechen-Backend.
//! - Mehrere Unbekannte: [`newton_system`], [`newton_system_with_jacobian`].
//! - Große lineare Systeme: [`conjugate_gradient`] (auch vorkonditioniert:
//!   [`preconditioned_conjugate_gradient`]), [`bicgstab()`], [`gmres()`].
//!
//! Alle Löser nehmen ihre Abbruchkriterien als [`Tolerance`] und liefern eine
//! [`Solution`]. Konvergiert ein Verfahren nicht, ist das ein
//! [`SolveError`], der den letzten Stand enthält.

mod bisection;
mod brent;
mod error;
mod iterative;
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
pub use iterative::{
    LinearResult, bicgstab, conjugate_gradient, gmres, gmres_work_len, preconditioned_conjugate_gradient,
};
pub use newton::{newton, newton_auto, newton_bracketed};
pub use secant::secant;
pub use solution::Solution;
pub use system::{SystemResult, newton_system, newton_system_with_jacobian};
pub use tolerance::Tolerance;
