//! Verfahren der verallgemeinerten minimalen Residuen (GMRES) mit Neustart.
//!
//! Aufbau nach Wikipedia, „Generalized minimal residual method“, Abschnitt
//! „Example code“: Arnoldi-Verfahren mit modifiziertem Gram-Schmidt und
//! Givens-Drehungen für das kleine Ausgleichsproblem.

use super::util::{LinearResult, axpy, dot, initial_residual, state};
use crate::linalg::LinearOperator;
use crate::scalar::RealElementary;
use crate::solve::{SolveError, Tolerance};

/// Nötige Länge des Arbeitsspeichers für [`gmres`] bei `size` Unbekannten
/// und Neustart nach `restart` Schritten.
///
/// `None`, wenn die Größe nicht darstellbar ist.
pub fn gmres_work_len(size: usize, restart: usize) -> Option<usize> {
    let columns = restart.checked_add(1)?;
    // Basis, Hilfsvektor, Hessenberg-Matrix, Drehungen und Lösung, rechte Seite.
    columns
        .checked_mul(size)?
        .checked_add(size)?
        .checked_add(columns.checked_mul(restart)?)?
        .checked_add(restart.checked_mul(3)?)?
        .checked_add(columns)
}

/// Spalte `index` der Basis.
#[inline]
fn column<T>(basis: &[T], size: usize, index: usize) -> &[T] {
    &basis[index * size..(index + 1) * size]
}

/// Löst `A·x = b` für eine allgemeine quadratische Matrix mit GMRES.
///
/// `x` enthält beim Aufruf den Startwert (z. B. null) und danach die Lösung.
/// Nach `restart` Schritten wird neu gestartet; größere Werte konvergieren
/// besser, brauchen aber mehr Speicher. `work` muss mindestens
/// [`gmres_work_len`]`(n, restart)` Elemente haben.
///
/// Fehler: [`SolveError::InvalidInput`] bei unpassenden Größen oder
/// `restart == 0`, [`SolveError::Singular`], wenn das Verfahren
/// zusammenbricht, sonst [`SolveError::NotConverged`] mit dem letzten
/// Restfehler.
// Die Schleifen brauchen den Index für mehrere Felder zugleich.
#[allow(clippy::needless_range_loop)]
pub fn gmres<T, A>(
    a: &A,
    b: &[T],
    x: &mut [T],
    restart: usize,
    work: &mut [T],
    tol: &Tolerance<T>,
) -> LinearResult<T>
where
    T: RealElementary,
    A: LinearOperator<T>,
{
    let size = b.len();
    if restart == 0 || a.rows() != size || a.cols() != size || x.len() != size {
        return Err(SolveError::InvalidInput);
    }
    let Some(needed) = gmres_work_len(size, restart) else {
        return Err(SolveError::InvalidInput);
    };
    let Some(work) = work.get_mut(..needed) else {
        return Err(SolveError::InvalidInput);
    };
    let stride = restart + 1;
    let (basis, rest) = work.split_at_mut(stride * size);
    let (w, rest) = rest.split_at_mut(size);
    // Hessenberg-Matrix spaltenweise: Element (zeile, spalte) liegt an
    // spalte · stride + zeile.
    let (h, rest) = rest.split_at_mut(stride * restart);
    let (cs, rest) = rest.split_at_mut(restart);
    let (sn, rest) = rest.split_at_mut(restart);
    let (y, g) = rest.split_at_mut(restart);

    let threshold = tol.abs_f + tol.rel_x * dot(b, b).sqrt();
    let mut iterations = 0;
    loop {
        initial_residual(a, b, x, w)?;
        let beta = dot(w, w).sqrt();
        if beta.is_nan() {
            return Err(SolveError::NonFinite(state(beta, iterations)));
        }
        if beta <= threshold {
            return Ok(state(beta, iterations));
        }
        if iterations >= tol.max_iterations {
            return Err(SolveError::NotConverged(state(beta, iterations)));
        }

        let inverse = beta.recip();
        for (target, &source) in basis[..size].iter_mut().zip(w.iter()) {
            *target = source * inverse;
        }
        g.fill(T::ZERO);
        g[0] = beta;

        let mut used = 0;
        for k in 0..restart {
            if iterations >= tol.max_iterations {
                break;
            }
            iterations += 1;

            // Arnoldi-Schritt: w = A·q_k, orthogonalisiert gegen die Basis.
            if a.apply(column(basis, size, k), w).is_err() {
                return Err(SolveError::InvalidInput);
            }
            let start = k * stride;
            for i in 0..=k {
                let q = column(basis, size, i);
                let coefficient = dot(w, q);
                h[start + i] = coefficient;
                axpy(-coefficient, q, w);
            }
            let next_norm = dot(w, w).sqrt();
            h[start + k + 1] = next_norm;
            if next_norm > T::ZERO {
                let inverse = next_norm.recip();
                let target = &mut basis[(k + 1) * size..(k + 2) * size];
                for (slot, &source) in target.iter_mut().zip(w.iter()) {
                    *slot = source * inverse;
                }
            }

            // Bisherige Givens-Drehungen auf die neue Spalte anwenden.
            for i in 0..k {
                let upper = cs[i] * h[start + i] + sn[i] * h[start + i + 1];
                h[start + i + 1] = cs[i] * h[start + i + 1] - sn[i] * h[start + i];
                h[start + i] = upper;
            }
            // Neue Drehung, die das Element unter der Diagonale auslöscht.
            let radius = h[start + k].hypot(h[start + k + 1]);
            if radius == T::ZERO || radius.is_nan() {
                return Err(SolveError::Singular(state(g[k].abs(), iterations)));
            }
            cs[k] = h[start + k] / radius;
            sn[k] = h[start + k + 1] / radius;
            h[start + k] = radius;
            h[start + k + 1] = T::ZERO;
            g[k + 1] = -sn[k] * g[k];
            g[k] = cs[k] * g[k];
            used = k + 1;

            // |g[k + 1]| ist der Restfehler nach diesem Schritt.
            if g[k + 1].abs() <= threshold || next_norm == T::ZERO {
                break;
            }
        }

        // Rückwärtseinsetzen im oberen Dreieck, dann x aktualisieren.
        for i in (0..used).rev() {
            let mut sum = g[i];
            for j in i + 1..used {
                sum -= h[j * stride + i] * y[j];
            }
            y[i] = sum / h[i * stride + i];
        }
        for j in 0..used {
            axpy(y[j], column(basis, size, j), x);
        }
    }
}
