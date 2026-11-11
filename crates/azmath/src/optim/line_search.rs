//! Gemeinsame Bausteine der gradientenbasierten Verfahren.

use crate::linalg::Vector;
use crate::scalar::RealElementary;
use crate::solve::{Solution, SolveError, SystemResult, Tolerance};

/// Höchstzahl der Halbierungen in der Liniensuche.
const MAX_HALVINGS: usize = 60;

/// Punkt mit Funktionswert und Gradient.
#[derive(Clone, Copy)]
pub(super) struct Iterate<T, const N: usize> {
    pub(super) x: Vector<T, N>,
    pub(super) value: T,
    pub(super) gradient: Vector<T, N>,
}

impl<T: RealElementary, const N: usize> Iterate<T, N> {
    pub(super) fn at(f: &mut impl FnMut(Vector<T, N>) -> (T, Vector<T, N>), x: Vector<T, N>) -> Self {
        let (value, gradient) = f(x);
        Self { x, value, gradient }
    }

    pub(super) fn state(&self, iterations: u32) -> Solution<Vector<T, N>, T> {
        Solution { x: self.x, residual: self.gradient.norm(), iterations }
    }
}

/// Rückwärts-Liniensuche nach Armijo entlang `direction`, beginnend mit der
/// Schrittweite `first`.
///
/// Liefert die angenommene Schrittweite und den neuen Punkt, oder `None`,
/// wenn sich der Funktionswert nicht ausreichend verkleinern ließ.
pub(super) fn backtrack<T, F, const N: usize>(
    f: &mut F,
    current: &Iterate<T, N>,
    direction: Vector<T, N>,
    first: T,
) -> Option<(T, Iterate<T, N>)>
where
    T: RealElementary,
    F: FnMut(Vector<T, N>) -> (T, Vector<T, N>),
{
    let slope = current.gradient.dot(direction);
    let sufficient = T::from_f64(1.0e-4);
    let mut step = first;
    for _ in 0..MAX_HALVINGS {
        let candidate = Iterate::at(f, current.x + direction * step);
        if candidate.value <= current.value + sufficient * step * slope {
            return Some((step, candidate));
        }
        step *= T::from_f64(0.5);
    }
    None
}

/// Prüft die Abbruchkriterien vor einem Schritt.
pub(super) fn check_start<T: RealElementary, const N: usize>(
    current: &Iterate<T, N>,
    iterations: u32,
    tol: &Tolerance<T>,
) -> Option<SystemResult<T, N>> {
    let state = current.state(iterations);
    if !current.value.is_finite() || !state.residual.is_finite() {
        return Some(Err(SolveError::NonFinite(state)));
    }
    if state.residual == T::ZERO || tol.f_converged(state.residual) {
        return Some(Ok(state));
    }
    None
}

/// Ergebnis, wenn die Liniensuche keinen Fortschritt mehr findet: Ist der
/// Gradient schon klein gegenüber dem Funktionswert, gilt die
/// Rechengenauigkeit als erreicht.
pub(super) fn stalled<T: RealElementary, const N: usize>(
    current: &Iterate<T, N>,
    iterations: u32,
) -> SystemResult<T, N> {
    let state = current.state(iterations);
    let working_precision = T::EPSILON.sqrt() * (T::ONE + current.value.abs());
    if state.residual <= working_precision {
        Ok(state)
    } else {
        Err(SolveError::NotConverged(state))
    }
}
