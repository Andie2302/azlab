//! Newton-Verfahren für nichtlineare Gleichungssysteme.

use super::error::SolveError;
use super::solution::Solution;
use super::tolerance::Tolerance;
use crate::autodiff::jacobian;
use crate::linalg::{Matrix, Vector};
use crate::number::DualN;
use crate::scalar::RealElementary;

/// Ergebnis der Löser für Gleichungssysteme mit `N` Unbekannten.
pub type SystemResult<T, const N: usize> =
    Result<Solution<Vector<T, N>, T>, SolveError<Vector<T, N>, T>>;

/// Höchstzahl der Halbierungen in der Liniensuche.
const MAX_HALVINGS: usize = 40;

/// Gedämpftes Newton-Verfahren für `F(x) = 0` mit `N` Gleichungen und `N`
/// Unbekannten; die Jacobi-Matrix liefert der Aufrufer.
///
/// `f` gibt das Paar `(F(x), Jacobi-Matrix)` zurück. Jeder Schritt löst ein
/// lineares System über die LU-Zerlegung. Verkleinert der volle Schritt den
/// Restfehler `‖F‖` nicht, wird er so lange halbiert, bis er es tut.
///
/// Lässt sich der Restfehler gar nicht mehr verkleinern, obwohl der
/// Newton-Schritt schon winzig ist, gilt die Rechengenauigkeit als erreicht
/// und das Ergebnis als konvergiert.
pub fn newton_system_with_jacobian<T, F, const N: usize>(
    mut f: F,
    x0: Vector<T, N>,
    tol: &Tolerance<T>,
) -> SystemResult<T, N>
where
    T: RealElementary,
    F: FnMut(Vector<T, N>) -> (Vector<T, N>, Matrix<T, N, N>),
{
    let half = T::from_f64(0.5);
    let mut x = x0;
    let (mut value, mut jacobi) = f(x);
    let mut residual = value.norm();

    for iteration in 1..=tol.max_iterations {
        let state = Solution { x, residual, iterations: iteration - 1 };
        if !residual.is_finite() {
            return Err(SolveError::NonFinite(state));
        }
        if residual == T::ZERO || tol.f_converged(residual) {
            return Ok(state);
        }
        let Ok(step) = jacobi.solve(value) else {
            return Err(SolveError::Singular(state));
        };

        // Liniensuche: Schritt halbieren, bis der Restfehler sinkt.
        let mut damping = T::ONE;
        let mut accepted = None;
        for _ in 0..MAX_HALVINGS {
            let candidate = x - step * damping;
            let (candidate_value, candidate_jacobi) = f(candidate);
            let candidate_residual = candidate_value.norm();
            if candidate_residual < residual {
                accepted = Some((candidate, candidate_value, candidate_jacobi, candidate_residual));
                break;
            }
            damping *= half;
        }
        let Some((next_x, next_value, next_jacobi, next_residual)) = accepted else {
            let working_precision = T::EPSILON.sqrt() * (T::ONE + x.norm());
            return if step.norm() <= working_precision {
                Ok(state)
            } else {
                Err(SolveError::NotConverged(state))
            };
        };

        let moved = (step * damping).norm();
        x = next_x;
        value = next_value;
        jacobi = next_jacobi;
        residual = next_residual;
        if tol.x_converged(moved, x.norm()) {
            return Ok(Solution { x, residual, iterations: iteration });
        }
    }
    Err(SolveError::NotConverged(Solution {
        x,
        residual,
        iterations: tol.max_iterations,
    }))
}

/// Gedämpftes Newton-Verfahren für `F(x) = 0` mit automatischer
/// Jacobi-Matrix.
///
/// `f` wird mit Dualzahlen aufgerufen; man schreibt es einmal generisch über
/// die Skalar-Traits. Sonst wie [`newton_system_with_jacobian`].
pub fn newton_system<T, F, const N: usize>(
    mut f: F,
    x0: Vector<T, N>,
    tol: &Tolerance<T>,
) -> SystemResult<T, N>
where
    T: RealElementary,
    F: FnMut(Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, N>,
{
    newton_system_with_jacobian(|x| jacobian(&mut f, x), x0, tol)
}
