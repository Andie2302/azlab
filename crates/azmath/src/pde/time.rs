//! Impliziter Zeitschritt für lineare Probleme `u' = A·u`.

use crate::linalg::{AffineOperator, LinearOperator};
use crate::scalar::RealElementary;
use crate::solve::{LinearResult, SolveError, Tolerance, bicgstab, conjugate_gradient};

/// Wahl des iterativen Lösers für das lineare System eines Zeitschritts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinearSolver {
    /// Konjugierte Gradienten: für symmetrische Operatoren wie den
    /// Laplace-Operator (Wärmeleitung, Diffusion).
    ConjugateGradient,
    /// BiCGSTAB: für unsymmetrische Operatoren (z. B. mit Transport).
    BiCgStab,
}

/// Ein Zeitschritt der Theta-Methode für `u' = A·u`:
/// `(I − theta·dt·A)·u_neu = (I + (1 − theta)·dt·A)·u_alt`.
///
/// `theta = 1` ist das implizite Euler-Verfahren (1. Ordnung, stark
/// dämpfend), `theta = 1/2` das Crank-Nicolson-Verfahren (2. Ordnung),
/// `theta = 0` das explizite Euler-Verfahren. `u` wird an Ort und Stelle
/// ersetzt. `work` braucht mindestens `7 · u.len()` Elemente.
///
/// Das Ergebnis ist das des linearen Lösers; `tol` steuert ihn.
pub fn theta_step<T, A>(
    a: &A,
    u: &mut [T],
    dt: T,
    theta: T,
    solver: LinearSolver,
    work: &mut [T],
    tol: &Tolerance<T>,
) -> LinearResult<T>
where
    T: RealElementary,
    A: LinearOperator<T>,
{
    let size = u.len();
    let Some(work) = size.checked_mul(7).and_then(|needed| work.get_mut(..needed)) else {
        return Err(SolveError::InvalidInput);
    };
    let (rhs, rest) = work.split_at_mut(size);
    let explicit = AffineOperator::new(a, T::ONE, (T::ONE - theta) * dt);
    if explicit.apply(u, rhs).is_err() {
        return Err(SolveError::InvalidInput);
    }
    let implicit = AffineOperator::new(a, T::ONE, -theta * dt);
    match solver {
        LinearSolver::ConjugateGradient => conjugate_gradient(&implicit, rhs, u, rest, tol),
        LinearSolver::BiCgStab => bicgstab(&implicit, rhs, u, rest, tol),
    }
}
