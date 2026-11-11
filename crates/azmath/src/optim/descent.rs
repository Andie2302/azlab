//! Gradientenabstieg.

use super::line_search::{Iterate, backtrack, check_start, stalled};
use crate::linalg::Vector;
use crate::scalar::RealElementary;
use crate::solve::{SolveError, SystemResult, Tolerance};

/// Minimiert `f` durch Gradientenabstieg mit Liniensuche ab `x0`.
///
/// `f` liefert `(Wert, Gradient)`. Jeder Schritt geht in Richtung des
/// steilsten Abstiegs; die Schrittweite wird so gewählt, dass der
/// Funktionswert ausreichend sinkt. Einfach und robust, bei schlecht
/// konditionierten Problemen aber langsam – dann [`bfgs`](super::bfgs).
pub fn gradient_descent<T, F, const N: usize>(
    mut f: F,
    x0: Vector<T, N>,
    tol: &Tolerance<T>,
) -> SystemResult<T, N>
where
    T: RealElementary,
    F: FnMut(Vector<T, N>) -> (T, Vector<T, N>),
{
    let mut current = Iterate::at(&mut f, x0);
    let mut step = T::ONE;
    for iteration in 1..=tol.max_iterations {
        if let Some(result) = check_start(&current, iteration - 1, tol) {
            return result;
        }
        let direction = -current.gradient;
        let Some((taken, next)) = backtrack(&mut f, &current, direction, step) else {
            return stalled(&current, iteration - 1);
        };
        let moved = (direction * taken).norm();
        current = next;
        // Nächster Versuch etwas mutiger als der letzte Erfolg.
        step = taken * T::from_f64(2.0);
        if tol.x_converged(moved, current.x.norm()) {
            return Ok(current.state(iteration));
        }
    }
    match check_start(&current, tol.max_iterations, tol) {
        Some(result) => result,
        None => Err(SolveError::NotConverged(current.state(tol.max_iterations))),
    }
}
