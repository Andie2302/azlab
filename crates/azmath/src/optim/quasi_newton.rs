//! Quasi-Newton-Verfahren: BFGS und L-BFGS.

use super::line_search::{Iterate, backtrack, check_start, stalled};
use crate::linalg::{Matrix, Vector};
use crate::scalar::RealElementary;
use crate::solve::{SolveError, SystemResult, Tolerance};

/// `true`, wenn das Paar aus Schritt `s` und Gradientenänderung `y` die
/// Krümmungsbedingung deutlich erfüllt; nur dann bleibt die Näherung der
/// Hesse-Matrix positiv definit.
fn has_curvature<T: RealElementary, const N: usize>(s: Vector<T, N>, y: Vector<T, N>) -> bool {
    s.dot(y) > T::EPSILON * s.norm() * y.norm()
}

/// Minimiert `f` mit dem BFGS-Verfahren ab `x0`.
///
/// `f` liefert `(Wert, Gradient)`. Aus den Gradienten aufeinanderfolgender
/// Schritte entsteht eine Näherung der inversen Hesse-Matrix; damit
/// konvergiert das Verfahren nahe dem Minimum sehr schnell. Speichert eine
/// `N × N`-Matrix; für sehr viele Unbekannte [`lbfgs`].
pub fn bfgs<T, F, const N: usize>(mut f: F, x0: Vector<T, N>, tol: &Tolerance<T>) -> SystemResult<T, N>
where
    T: RealElementary,
    F: FnMut(Vector<T, N>) -> (T, Vector<T, N>),
{
    let identity = Matrix::<T, N, N>::identity();
    let mut inverse_hessian = identity;
    let mut current = Iterate::at(&mut f, x0);
    for iteration in 1..=tol.max_iterations {
        if let Some(result) = check_start(&current, iteration - 1, tol) {
            return result;
        }
        let mut direction = -(inverse_hessian * current.gradient);
        if current.gradient.dot(direction) >= T::ZERO {
            // Keine Abstiegsrichtung mehr: Näherung verwerfen.
            inverse_hessian = identity;
            direction = -current.gradient;
        }
        let Some((taken, next)) = backtrack(&mut f, &current, direction, T::ONE) else {
            if inverse_hessian == identity {
                return stalled(&current, iteration - 1);
            }
            inverse_hessian = identity;
            continue;
        };
        let s = direction * taken;
        let y = next.gradient - current.gradient;
        if has_curvature(s, y) {
            // H ← (I − ρ·s·yᵀ)·H·(I − ρ·y·sᵀ) + ρ·s·sᵀ
            let rho = s.dot(y).recip();
            let left = identity - s.outer(y) * rho;
            inverse_hessian = left * inverse_hessian * left.transpose() + s.outer(s) * rho;
        }
        current = next;
        if tol.x_converged(s.norm(), current.x.norm()) {
            return Ok(current.state(iteration));
        }
    }
    match check_start(&current, tol.max_iterations, tol) {
        Some(result) => result,
        None => Err(SolveError::NotConverged(current.state(tol.max_iterations))),
    }
}

/// Minimiert `f` mit dem L-BFGS-Verfahren ab `x0`.
///
/// Wie [`bfgs`], merkt sich aber statt einer Matrix nur die letzten `M`
/// Schritte (`M >= 1`, üblich sind 5 bis 20). Der Speicherbedarf wächst
/// damit nur linear mit der Zahl der Unbekannten.
pub fn lbfgs<T, F, const N: usize, const M: usize>(
    mut f: F,
    x0: Vector<T, N>,
    tol: &Tolerance<T>,
) -> SystemResult<T, N>
where
    T: RealElementary,
    F: FnMut(Vector<T, N>) -> (T, Vector<T, N>),
{
    if M == 0 {
        return Err(SolveError::InvalidInput);
    }
    // Ringpuffer der letzten Paare (s, y) mit ρ = 1/(s·y).
    let mut steps = [Vector::<T, N>::ZERO; M];
    let mut changes = [Vector::<T, N>::ZERO; M];
    let mut rhos = [T::ZERO; M];
    let mut stored = 0;
    let mut newest = 0;

    let mut current = Iterate::at(&mut f, x0);
    for iteration in 1..=tol.max_iterations {
        if let Some(result) = check_start(&current, iteration - 1, tol) {
            return result;
        }
        // Zwei-Schleifen-Rekursion: wendet die Näherung der inversen
        // Hesse-Matrix auf den Gradienten an.
        let mut q = current.gradient;
        let mut alphas = [T::ZERO; M];
        for back in 0..stored {
            let index = (newest + M - back) % M;
            alphas[index] = rhos[index] * steps[index].dot(q);
            q -= changes[index] * alphas[index];
        }
        if stored > 0 {
            let y = changes[newest];
            q *= steps[newest].dot(y) / y.dot(y);
        }
        for back in (0..stored).rev() {
            let index = (newest + M - back) % M;
            let beta = rhos[index] * changes[index].dot(q);
            q += steps[index] * (alphas[index] - beta);
        }
        let mut direction = -q;
        if current.gradient.dot(direction) >= T::ZERO {
            stored = 0;
            direction = -current.gradient;
        }

        let Some((taken, next)) = backtrack(&mut f, &current, direction, T::ONE) else {
            if stored == 0 {
                return stalled(&current, iteration - 1);
            }
            stored = 0;
            continue;
        };
        let s = direction * taken;
        let y = next.gradient - current.gradient;
        if has_curvature(s, y) {
            newest = if stored == 0 { 0 } else { (newest + 1) % M };
            steps[newest] = s;
            changes[newest] = y;
            rhos[newest] = s.dot(y).recip();
            stored = (stored + 1).min(M);
        }
        current = next;
        if tol.x_converged(s.norm(), current.x.norm()) {
            return Ok(current.state(iteration));
        }
    }
    match check_start(&current, tol.max_iterations, tol) {
        Some(result) => result,
        None => Err(SolveError::NotConverged(current.state(tol.max_iterations))),
    }
}
