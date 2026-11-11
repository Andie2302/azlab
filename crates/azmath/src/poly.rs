//! Polynome: Auswertung und Nullstellen.
//!
//! Koeffizienten stehen immer aufsteigend: `c[0] + c[1]·x + c[2]·x² + …`.
//!
//! | Aufgabe | Werkzeug |
//! |---------|----------|
//! | Auswerten | [`Polynomial::eval`], [`horner`], [`clenshaw`] |
//! | Nullstellen bis Grad 4 | [`solve_linear`] bis [`solve_quartic`] |
//! | Alle komplexen Nullstellen, jeder Grad | [`complex_roots`] (ohne Garantie) |
//! | Reelle Nullstellen in einem Intervall, lückenlos | [`isolate_real_roots`], [`refine_root`] |

mod aberth;
mod closed;
mod eval;
mod isolate;
mod ops;
mod real_roots;
mod repr;

#[cfg(test)]
mod tests;

pub use aberth::complex_roots;
pub use closed::{
    solve_cubic, solve_linear, solve_quadratic, solve_quadratic_complex, solve_quartic,
};
pub use eval::{clenshaw, horner, horner_with_derivative};
pub use isolate::{BracketKind, RootBracket, isolate_real_roots, refine_root};
pub use real_roots::RealRoots;
pub use repr::Polynomial;
