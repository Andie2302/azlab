//! Stabilisiertes Verfahren der bikonjugierten Gradienten (BiCGSTAB).

use super::util::{LinearResult, axpy, dot, initial_residual, state};
use crate::linalg::LinearOperator;
use crate::scalar::RealElementary;
use crate::solve::{SolveError, Tolerance};

/// Löst `A·x = b` für eine allgemeine quadratische Matrix mit BiCGSTAB.
///
/// `x` enthält beim Aufruf den Startwert (z. B. null) und danach die Lösung.
/// `work` muss mindestens `6·n` Elemente haben.
///
/// Fehler: [`SolveError::InvalidInput`] bei unpassenden Größen,
/// [`SolveError::Singular`], wenn das Verfahren zusammenbricht, sonst
/// [`SolveError::NotConverged`] mit dem letzten Restfehler.
pub fn bicgstab<T, A>(
    a: &A,
    b: &[T],
    x: &mut [T],
    work: &mut [T],
    tol: &Tolerance<T>,
) -> LinearResult<T>
where
    T: RealElementary,
    A: LinearOperator<T>,
{
    let size = b.len();
    let Some(needed) = size.checked_mul(6) else {
        return Err(SolveError::InvalidInput);
    };
    let Some(work) = work.get_mut(..needed) else {
        return Err(SolveError::InvalidInput);
    };
    let (r, rest) = work.split_at_mut(size);
    let (shadow, rest) = rest.split_at_mut(size);
    let (p, rest) = rest.split_at_mut(size);
    let (v, rest) = rest.split_at_mut(size);
    let (s, t) = rest.split_at_mut(size);

    initial_residual(a, b, x, r)?;
    shadow.copy_from_slice(r);
    p.fill(T::ZERO);
    v.fill(T::ZERO);
    let threshold = tol.abs_f + tol.rel_x * dot(b, b).sqrt();
    let mut rho = T::ONE;
    let mut alpha = T::ONE;
    let mut omega = T::ONE;
    let mut residual = dot(r, r).sqrt();

    for iteration in 1..=tol.max_iterations {
        if residual.is_nan() {
            return Err(SolveError::NonFinite(state(residual, iteration - 1)));
        }
        if residual <= threshold {
            return Ok(state(residual, iteration - 1));
        }
        let breakdown = SolveError::Singular(state(residual, iteration - 1));

        let next_rho = dot(shadow, r);
        if next_rho == T::ZERO || omega == T::ZERO {
            return Err(breakdown);
        }
        let beta = (next_rho / rho) * (alpha / omega);
        // p = r + beta·(p − omega·v)
        for ((direction, &rest), &image) in p.iter_mut().zip(r.iter()).zip(v.iter()) {
            *direction = rest + beta * (*direction - omega * image);
        }
        if a.apply(p, v).is_err() {
            return Err(SolveError::InvalidInput);
        }
        let denominator = dot(shadow, v);
        if denominator == T::ZERO {
            return Err(breakdown);
        }
        alpha = next_rho / denominator;

        // s = r − alpha·v
        s.copy_from_slice(r);
        axpy(-alpha, v, s);
        let half_step = dot(s, s).sqrt();
        if half_step <= threshold {
            axpy(alpha, p, x);
            return Ok(state(half_step, iteration));
        }

        if a.apply(s, t).is_err() {
            return Err(SolveError::InvalidInput);
        }
        let tt = dot(t, t);
        if tt == T::ZERO {
            return Err(breakdown);
        }
        omega = dot(t, s) / tt;

        axpy(alpha, p, x);
        axpy(omega, s, x);
        // r = s − omega·t
        r.copy_from_slice(s);
        axpy(-omega, t, r);
        rho = next_rho;
        residual = dot(r, r).sqrt();
    }

    if residual <= threshold {
        Ok(state(residual, tol.max_iterations))
    } else {
        Err(SolveError::NotConverged(state(residual, tol.max_iterations)))
    }
}
