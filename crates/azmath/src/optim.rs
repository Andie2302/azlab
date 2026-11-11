//! Optimierung: Minima von Funktionen mehrerer Veränderlicher finden.
//!
//! | Verfahren | Braucht | Wofür |
//! |-----------|---------|-------|
//! | [`gradient_descent`] | Gradient | einfach und robust, aber langsam |
//! | [`bfgs`] | Gradient | Standard für glatte Probleme mäßiger Größe |
//! | [`lbfgs`] | Gradient | wie BFGS, mit wenig Speicher für viele Unbekannte |
//! | [`nelder_mead`] | nur Funktionswerte | nicht glatte oder verrauschte Funktionen |
//! | [`levenberg_marquardt`] | Residuen und Jacobi-Matrix | Ausgleichsrechnung (Kurvenanpassung) |
//!
//! Die gradientenbasierten Verfahren erwarten eine Funktion, die
//! `(Wert, Gradient)` liefert. Den Gradienten bekommt man ohne Handarbeit
//! aus [`autodiff::gradient`](crate::autodiff::gradient) oder
//! [`autodiff::reverse_gradient`](crate::autodiff::reverse_gradient):
//! `bfgs(|x| gradient(|p| f(p), x), x0, &tol)`.
//!
//! Abbruch: Die Norm des Gradienten ist höchstens `abs_f`, oder der letzte
//! Schritt ist nach [`Tolerance`](crate::solve::Tolerance) klein genug. Im
//! Ergebnis steht unter `residual` die Norm des Gradienten (bei
//! [`nelder_mead`] die Spanne der Funktionswerte, bei
//! [`levenberg_marquardt`] die Norm der Residuen).

mod descent;
mod least_squares;
mod line_search;
mod quasi_newton;
mod simplex;

#[cfg(test)]
mod tests;

pub use descent::gradient_descent;
pub use least_squares::{levenberg_marquardt, levenberg_marquardt_auto};
pub use quasi_newton::{bfgs, lbfgs};
pub use simplex::{nelder_mead, simplex_around};
