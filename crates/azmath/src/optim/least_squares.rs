//! Nichtlineare Ausgleichsrechnung mit dem Levenberg-Marquardt-Verfahren.

use crate::autodiff::jacobian;
use crate::linalg::{Matrix, Vector};
use crate::number::DualN;
use crate::scalar::RealElementary;
use crate::solve::{Solution, SolveError, SystemResult, Tolerance};

/// Minimiert die Summe der Quadrate der `M` Residuen `r(x)` über `N`
/// Unbekannte mit dem Levenberg-Marquardt-Verfahren.
///
/// `f` liefert `(Residuen, Jacobi-Matrix)`. Das Verfahren wechselt
/// stufenlos zwischen Gauß-Newton (schnell nahe der Lösung) und
/// Gradientenabstieg (sicher weit weg davon). Typische Anwendung:
/// Parameter eines Modells an Messwerte anpassen; dann ist jedes Residuum
/// „Modell minus Messwert“.
///
/// Abbruch, wenn die Norm des Gradienten `Jᵀ·r` höchstens `abs_f` ist oder
/// der Schritt klein genug. Im Ergebnis ist `residual` die Norm der
/// Residuen.
pub fn levenberg_marquardt<T, F, const N: usize, const M: usize>(
    mut f: F,
    x0: Vector<T, N>,
    tol: &Tolerance<T>,
) -> SystemResult<T, N>
where
    T: RealElementary,
    F: FnMut(Vector<T, N>) -> (Vector<T, M>, Matrix<T, M, N>),
{
    let ten = T::from_f64(10.0);
    let (smallest, largest) = (T::from_f64(1.0e-12), T::from_f64(1.0e12));
    let mut damping = T::from_f64(1.0e-3);
    let mut x = x0;
    let (mut residuals, mut jacobi) = f(x);
    let mut cost = residuals.norm_sqr();

    for iteration in 1..=tol.max_iterations {
        let state = Solution { x, residual: cost.sqrt(), iterations: iteration - 1 };
        if !cost.is_finite() {
            return Err(SolveError::NonFinite(state));
        }
        let transposed = jacobi.transpose();
        let gradient = transposed * residuals;
        let gradient_norm = gradient.norm();
        if gradient_norm == T::ZERO || tol.f_converged(gradient_norm) {
            return Ok(state);
        }
        // Gedämpfte Normalgleichungen (JᵀJ + λ·D)·δ = −Jᵀr; D ist die
        // Diagonale von JᵀJ, damit die Dämpfung zur Skalierung passt.
        let normal = transposed * jacobi;
        let mut damped = normal;
        for index in 0..N {
            damped[(index, index)] += damping * normal[(index, index)].max(T::EPSILON);
        }
        let accepted = damped.solve(-gradient).ok().and_then(|step| {
            let candidate = x + step;
            let (next_residuals, next_jacobi) = f(candidate);
            let next_cost = next_residuals.norm_sqr();
            (next_cost < cost).then_some((step, candidate, next_residuals, next_jacobi, next_cost))
        });
        match accepted {
            Some((step, candidate, next_residuals, next_jacobi, next_cost)) => {
                x = candidate;
                residuals = next_residuals;
                jacobi = next_jacobi;
                cost = next_cost;
                damping = (damping / ten).max(smallest);
                if tol.x_converged(step.norm(), x.norm()) {
                    return Ok(Solution { x, residual: cost.sqrt(), iterations: iteration });
                }
            }
            None => {
                damping *= ten;
                if damping > largest {
                    // Kein Fortschritt mehr möglich: Rechengenauigkeit
                    // erreicht, wenn der Gradient klein gegen die Kosten ist.
                    let working_precision = T::EPSILON.sqrt() * (T::ONE + cost);
                    return if gradient_norm <= working_precision {
                        Ok(state)
                    } else {
                        Err(SolveError::NotConverged(state))
                    };
                }
            }
        }
    }
    Err(SolveError::NotConverged(Solution {
        x,
        residual: cost.sqrt(),
        iterations: tol.max_iterations,
    }))
}

/// Wie [`levenberg_marquardt`], mit automatischer Jacobi-Matrix.
///
/// `f` wird mit Dualzahlen aufgerufen und liefert nur die Residuen.
pub fn levenberg_marquardt_auto<T, F, const N: usize, const M: usize>(
    mut f: F,
    x0: Vector<T, N>,
    tol: &Tolerance<T>,
) -> SystemResult<T, N>
where
    T: RealElementary,
    F: FnMut(Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, M>,
{
    levenberg_marquardt(|x| jacobian(&mut f, x), x0, tol)
}
