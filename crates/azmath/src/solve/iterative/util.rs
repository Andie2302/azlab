//! Gemeinsame Bausteine der iterativen Löser.

use crate::linalg::LinearOperator;
use crate::scalar::Real;
use crate::solve::{Solution, SolveError};

/// Ergebnis der iterativen Löser: Die Lösung steht im Vektor `x` des
/// Aufrufers, zurück kommen Restfehler und Iterationszahl.
pub type LinearResult<T> = Result<Solution<(), T>, SolveError<(), T>>;

/// Skalarprodukt zweier gleich langer Slices.
#[inline]
pub(super) fn dot<T: Real>(a: &[T], b: &[T]) -> T {
    a.iter().zip(b).fold(T::ZERO, |acc, (&x, &y)| acc + x * y)
}

/// `y += alpha · x`.
#[inline]
pub(super) fn axpy<T: Real>(alpha: T, x: &[T], y: &mut [T]) {
    for (target, &source) in y.iter_mut().zip(x) {
        *target += alpha * source;
    }
}

/// Zustand für Fehler- und Erfolgsmeldungen.
#[inline]
pub(super) fn state<T>(residual: T, iterations: u32) -> Solution<(), T> {
    Solution { x: (), residual, iterations }
}

/// Prüft die Größen und berechnet den Anfangsrest `r = b − A·x`.
pub(super) fn initial_residual<T, A>(
    a: &A,
    b: &[T],
    x: &[T],
    r: &mut [T],
) -> Result<(), SolveError<(), T>>
where
    T: Real,
    A: LinearOperator<T>,
{
    let size = b.len();
    if a.rows() != size || a.cols() != size || x.len() != size || r.len() != size {
        return Err(SolveError::InvalidInput);
    }
    if a.apply(x, r).is_err() {
        return Err(SolveError::InvalidInput);
    }
    for (rest, &rhs) in r.iter_mut().zip(b) {
        *rest = rhs - *rest;
    }
    Ok(())
}
