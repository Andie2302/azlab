//! Verfahren der konjugierten Gradienten.

use super::util::{LinearResult, axpy, dot, initial_residual, state};
use crate::linalg::LinearOperator;
use crate::scalar::RealElementary;
use crate::solve::{SolveError, Tolerance};

/// Löst `A·x = b` für eine symmetrisch positiv definite Matrix mit dem
/// Verfahren der konjugierten Gradienten.
///
/// `x` enthält beim Aufruf den Startwert (z. B. null) und danach die Lösung.
/// `work` muss mindestens `3·n` Elemente haben.
///
/// Fehler: [`SolveError::InvalidInput`] bei unpassenden Größen,
/// [`SolveError::Singular`], wenn die Matrix nicht positiv definit ist,
/// sonst [`SolveError::NotConverged`] mit dem letzten Restfehler.
pub fn conjugate_gradient<T, A>(
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
    let Some(needed) = size.checked_mul(3) else {
        return Err(SolveError::InvalidInput);
    };
    let Some(work) = work.get_mut(..needed) else {
        return Err(SolveError::InvalidInput);
    };
    let (r, rest) = work.split_at_mut(size);
    let (p, ap) = rest.split_at_mut(size);

    initial_residual(a, b, x, r)?;
    p.copy_from_slice(r);
    let threshold = tol.abs_f + tol.rel_x * dot(b, b).sqrt();
    let mut rr = dot(r, r);

    for iteration in 1..=tol.max_iterations {
        let residual = rr.sqrt();
        if residual.is_nan() {
            return Err(SolveError::NonFinite(state(residual, iteration - 1)));
        }
        if residual <= threshold {
            return Ok(state(residual, iteration - 1));
        }
        if a.apply(p, ap).is_err() {
            return Err(SolveError::InvalidInput);
        }
        let curvature = dot(p, ap);
        if curvature <= T::ZERO || curvature.is_nan() {
            return Err(SolveError::Singular(state(residual, iteration - 1)));
        }
        let alpha = rr / curvature;
        axpy(alpha, p, x);
        axpy(-alpha, ap, r);
        let next_rr = dot(r, r);
        let beta = next_rr / rr;
        for (direction, &rest) in p.iter_mut().zip(r.iter()) {
            *direction = rest + beta * *direction;
        }
        rr = next_rr;
    }

    let residual = rr.sqrt();
    if residual <= threshold {
        Ok(state(residual, tol.max_iterations))
    } else {
        Err(SolveError::NotConverged(state(residual, tol.max_iterations)))
    }
}

/// Wie [`conjugate_gradient`], aber mit Vorkonditionierer: `preconditioner`
/// wendet eine Näherung `M⁻¹ ≈ A⁻¹` an, die selbst symmetrisch positiv
/// definit sein muss (z. B. die Kehrwerte der Diagonale oder eine
/// unvollständige Cholesky-Zerlegung).
///
/// Ein guter Vorkonditionierer senkt die Zahl der Iterationen deutlich; für
/// Gitterprobleme wächst sie sonst mit der Auflösung. Das Abbruchkriterium
/// misst weiterhin den echten Restfehler `‖b − A·x‖`.
///
/// `work` muss mindestens `4·n` Elemente haben.
///
/// Fehler wie bei [`conjugate_gradient`]; [`SolveError::Singular`] auch dann,
/// wenn der Vorkonditionierer nicht positiv definit ist.
pub fn preconditioned_conjugate_gradient<T, A, M>(
    a: &A,
    preconditioner: &M,
    b: &[T],
    x: &mut [T],
    work: &mut [T],
    tol: &Tolerance<T>,
) -> LinearResult<T>
where
    T: RealElementary,
    A: LinearOperator<T>,
    M: LinearOperator<T>,
{
    let size = b.len();
    let Some(needed) = size.checked_mul(4) else {
        return Err(SolveError::InvalidInput);
    };
    let Some(work) = work.get_mut(..needed) else {
        return Err(SolveError::InvalidInput);
    };
    let (r, rest) = work.split_at_mut(size);
    let (z, rest) = rest.split_at_mut(size);
    let (p, ap) = rest.split_at_mut(size);

    initial_residual(a, b, x, r)?;
    if preconditioner.rows() != size || preconditioner.cols() != size || preconditioner.apply(r, z).is_err() {
        return Err(SolveError::InvalidInput);
    }
    p.copy_from_slice(z);
    let threshold = tol.abs_f + tol.rel_x * dot(b, b).sqrt();
    let mut rz = dot(r, z);
    let mut residual = dot(r, r).sqrt();

    for iteration in 1..=tol.max_iterations {
        if residual.is_nan() {
            return Err(SolveError::NonFinite(state(residual, iteration - 1)));
        }
        if residual <= threshold {
            return Ok(state(residual, iteration - 1));
        }
        if a.apply(p, ap).is_err() {
            return Err(SolveError::InvalidInput);
        }
        let curvature = dot(p, ap);
        if curvature <= T::ZERO || curvature.is_nan() || rz <= T::ZERO || rz.is_nan() {
            return Err(SolveError::Singular(state(residual, iteration - 1)));
        }
        let alpha = rz / curvature;
        axpy(alpha, p, x);
        axpy(-alpha, ap, r);
        residual = dot(r, r).sqrt();
        if preconditioner.apply(r, z).is_err() {
            return Err(SolveError::InvalidInput);
        }
        let next_rz = dot(r, z);
        let beta = next_rz / rz;
        for (direction, &step) in p.iter_mut().zip(z.iter()) {
            *direction = step + beta * *direction;
        }
        rz = next_rz;
    }

    if residual <= threshold {
        Ok(state(residual, tol.max_iterations))
    } else {
        Err(SolveError::NotConverged(state(residual, tol.max_iterations)))
    }
}
