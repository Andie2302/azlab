//! Geometrisches Mehrgitterverfahren für die Poisson-Gleichung.
//!
//! Gelöst wird `−Δu = f` mit Wert null auf dem Rand, auf Gittern mit
//! `2^k − 1` inneren Punkten je Richtung. Ein V-Zyklus glättet mit
//! Gauß-Seidel, überträgt den Rest auf das gröbere Gitter (volle
//! Gewichtung), löst dort rekursiv und interpoliert die Korrektur zurück.
//! Der Aufwand je Zyklus ist proportional zur Zahl der Unbekannten, die
//! Zahl der Zyklen hängt nicht von der Gittergröße ab.

use crate::scalar::RealElementary;
use crate::solve::{LinearResult, Solution, SolveError, Tolerance};

/// Glättungsschritte vor und nach der Grobgitterkorrektur.
const SMOOTHING: usize = 2;

/// `true`, wenn `points + 1` eine Zweierpotenz ist (also `points = 2^k − 1`).
fn valid_size(points: usize) -> bool {
    points >= 1 && (points + 1).is_power_of_two()
}

fn norm<T: RealElementary>(values: &[T]) -> T {
    values.iter().fold(T::ZERO, |sum, &v| sum + v * v).sqrt()
}

/// Wert an Stelle `index`, null außerhalb (Randbedingung).
#[inline]
fn at<T: RealElementary>(values: &[T], index: Option<usize>) -> T {
    index.and_then(|i| values.get(i)).copied().unwrap_or(T::ZERO)
}

// ---------------------------------------------------------------- 1D ----

fn smooth_1d<T: RealElementary>(u: &mut [T], f: &[T], h2: T) {
    let half = T::from_f64(0.5);
    for i in 0..u.len() {
        let neighbours = at(u, i.checked_sub(1)) + at(u, Some(i + 1));
        u[i] = (f[i] * h2 + neighbours) * half;
    }
}

fn residual_1d<T: RealElementary>(u: &[T], f: &[T], h2: T, r: &mut [T]) {
    let two = T::from_f64(2.0);
    for (i, out) in r.iter_mut().enumerate() {
        let neighbours = at(u, i.checked_sub(1)) + at(u, Some(i + 1));
        *out = f[i] - (two * u[i] - neighbours) / h2;
    }
}

fn cycle_1d<T: RealElementary>(u: &mut [T], f: &[T], h2: T, work: &mut [T]) {
    let fine = u.len();
    if fine == 1 {
        u[0] = f[0] * h2 * T::from_f64(0.5);
        return;
    }
    let coarse = (fine - 1) / 2;
    let (r, rest) = work.split_at_mut(fine);
    let (coarse_f, rest) = rest.split_at_mut(coarse);
    let (coarse_u, rest) = rest.split_at_mut(coarse);

    for _ in 0..SMOOTHING {
        smooth_1d(u, f, h2);
    }
    residual_1d(u, f, h2, r);
    let (quarter, half) = (T::from_f64(0.25), T::from_f64(0.5));
    for (j, out) in coarse_f.iter_mut().enumerate() {
        *out = (r[2 * j] + r[2 * j + 2]) * quarter + r[2 * j + 1] * half;
    }
    coarse_u.fill(T::ZERO);
    cycle_1d(coarse_u, coarse_f, h2 * T::from_f64(4.0), rest);
    for i in 0..fine {
        // Ungerade Feinpunkte liegen auf Grobpunkten, gerade dazwischen.
        let correction = if i % 2 == 1 {
            coarse_u[i / 2]
        } else {
            (at(coarse_u, (i / 2).checked_sub(1)) + at(coarse_u, Some(i / 2))) * half
        };
        u[i] += correction;
    }
    for _ in 0..SMOOTHING {
        smooth_1d(u, f, h2);
    }
}

/// Löst `−u'' = f` mit `u = 0` an beiden Rändern auf `f.len() = 2^k − 1`
/// inneren Punkten im Abstand `h`.
///
/// `u` enthält den Startwert und danach die Lösung. `work` braucht
/// mindestens `5 · f.len()` Elemente. Konvergiert, wenn die Norm des Rests
/// höchstens `abs_f + rel_x · ‖f‖` ist; `max_iterations` begrenzt die Zahl
/// der V-Zyklen.
pub fn poisson_multigrid_1d<T: RealElementary>(
    f: &[T],
    u: &mut [T],
    h: T,
    work: &mut [T],
    tol: &Tolerance<T>,
) -> LinearResult<T> {
    let size = f.len();
    if !valid_size(size) || u.len() != size {
        return Err(SolveError::InvalidInput);
    }
    let Some(work) = size.checked_mul(5).and_then(|needed| work.get_mut(..needed)) else {
        return Err(SolveError::InvalidInput);
    };
    let h2 = h * h;
    let threshold = tol.abs_f + tol.rel_x * norm(f);
    let (r, rest) = work.split_at_mut(size);
    let mut iterations = 0;
    loop {
        residual_1d(u, f, h2, r);
        let residual = norm(r);
        let state = Solution { x: (), residual, iterations };
        if residual.is_nan() {
            return Err(SolveError::NonFinite(state));
        }
        if residual <= threshold {
            return Ok(state);
        }
        if iterations >= tol.max_iterations {
            return Err(SolveError::NotConverged(state));
        }
        cycle_1d(u, f, h2, rest);
        iterations += 1;
    }
}

// ---------------------------------------------------------------- 2D ----

/// Wert am Gitterpunkt `(row, col)` eines `n × n`-Felds, null außerhalb.
#[inline]
fn at_2d<T: RealElementary>(values: &[T], n: usize, row: Option<usize>, col: Option<usize>) -> T {
    match (row, col) {
        (Some(row), Some(col)) if row < n && col < n => values[row * n + col],
        _ => T::ZERO,
    }
}

fn cross_2d<T: RealElementary>(u: &[T], n: usize, row: usize, col: usize) -> T {
    at_2d(u, n, row.checked_sub(1), Some(col))
        + at_2d(u, n, Some(row + 1), Some(col))
        + at_2d(u, n, Some(row), col.checked_sub(1))
        + at_2d(u, n, Some(row), Some(col + 1))
}

fn smooth_2d<T: RealElementary>(u: &mut [T], f: &[T], n: usize, h2: T) {
    let quarter = T::from_f64(0.25);
    for row in 0..n {
        for col in 0..n {
            u[row * n + col] = (f[row * n + col] * h2 + cross_2d(u, n, row, col)) * quarter;
        }
    }
}

fn residual_2d<T: RealElementary>(u: &[T], f: &[T], n: usize, h2: T, r: &mut [T]) {
    let four = T::from_f64(4.0);
    for row in 0..n {
        for col in 0..n {
            let index = row * n + col;
            r[index] = f[index] - (four * u[index] - cross_2d(u, n, row, col)) / h2;
        }
    }
}

fn cycle_2d<T: RealElementary>(u: &mut [T], f: &[T], n: usize, h2: T, work: &mut [T]) {
    if n == 1 {
        u[0] = f[0] * h2 * T::from_f64(0.25);
        return;
    }
    let coarse = (n - 1) / 2;
    let (r, rest) = work.split_at_mut(n * n);
    let (coarse_f, rest) = rest.split_at_mut(coarse * coarse);
    let (coarse_u, rest) = rest.split_at_mut(coarse * coarse);

    for _ in 0..SMOOTHING {
        smooth_2d(u, f, n, h2);
    }
    residual_2d(u, f, n, h2, r);
    // Volle Gewichtung: Mitte 1/4, Kanten 1/8, Ecken 1/16.
    let (centre, edge, corner) = (T::from_f64(0.25), T::from_f64(0.125), T::from_f64(0.0625));
    for row in 0..coarse {
        for col in 0..coarse {
            let (i, j) = (2 * row + 1, 2 * col + 1);
            let value = |di: usize, dj: usize| r[(i + di - 1) * n + (j + dj - 1)];
            coarse_f[row * coarse + col] = value(1, 1) * centre
                + (value(0, 1) + value(2, 1) + value(1, 0) + value(1, 2)) * edge
                + (value(0, 0) + value(0, 2) + value(2, 0) + value(2, 2)) * corner;
        }
    }
    coarse_u.fill(T::ZERO);
    cycle_2d(coarse_u, coarse_f, coarse, h2 * T::from_f64(4.0), rest);

    // Bilineare Interpolation der Korrektur.
    let half = T::from_f64(0.5);
    let along = |index: usize| -> (Option<usize>, Option<usize>, T) {
        // Grobpunkte, zwischen denen der Feinpunkt liegt, und ihr Gewicht.
        if index % 2 == 1 {
            (Some(index / 2), None, T::ONE)
        } else {
            ((index / 2).checked_sub(1), Some(index / 2), half)
        }
    };
    for row in 0..n {
        let (row_a, row_b, row_weight) = along(row);
        for col in 0..n {
            let (col_a, col_b, col_weight) = along(col);
            let mut correction = T::ZERO;
            for r_index in [row_a, row_b] {
                for c_index in [col_a, col_b] {
                    correction += at_2d(coarse_u, coarse, r_index, c_index);
                }
            }
            u[row * n + col] += correction * row_weight * col_weight;
        }
    }
    for _ in 0..SMOOTHING {
        smooth_2d(u, f, n, h2);
    }
}

/// Löst `−Δu = f` mit `u = 0` auf dem Rand auf einem quadratischen Gitter
/// mit `n × n` inneren Punkten (`n = 2^k − 1`) im Abstand `h`.
///
/// `f` und `u` haben `n²` Elemente, zeilenweise. `u` enthält den Startwert
/// und danach die Lösung. `work` braucht mindestens `5 · n²` Elemente.
/// Abbruchkriterium wie bei [`poisson_multigrid_1d`].
pub fn poisson_multigrid_2d<T: RealElementary>(
    f: &[T],
    u: &mut [T],
    n: usize,
    h: T,
    work: &mut [T],
    tol: &Tolerance<T>,
) -> LinearResult<T> {
    let Some(size) = n.checked_mul(n) else {
        return Err(SolveError::InvalidInput);
    };
    if !valid_size(n) || f.len() != size || u.len() != size {
        return Err(SolveError::InvalidInput);
    }
    let Some(work) = size.checked_mul(5).and_then(|needed| work.get_mut(..needed)) else {
        return Err(SolveError::InvalidInput);
    };
    let h2 = h * h;
    let threshold = tol.abs_f + tol.rel_x * norm(f);
    let (r, rest) = work.split_at_mut(size);
    let mut iterations = 0;
    loop {
        residual_2d(u, f, n, h2, r);
        let residual = norm(r);
        let state = Solution { x: (), residual, iterations };
        if residual.is_nan() {
            return Err(SolveError::NonFinite(state));
        }
        if residual <= threshold {
            return Ok(state);
        }
        if iterations >= tol.max_iterations {
            return Err(SolveError::NotConverged(state));
        }
        cycle_2d(u, f, n, h2, rest);
        iterations += 1;
    }
}
