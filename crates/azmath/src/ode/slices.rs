//! Einzelschritte auf Slices für Systeme mit Größe zur Laufzeit.

use crate::linalg::LinalgError;
use crate::scalar::Real;

/// Ein Schritt des expliziten Euler-Verfahrens an Ort und Stelle.
///
/// `f(t, y, dy)` schreibt die Ableitung nach `dy`. `work` braucht mindestens
/// `y.len()` Elemente.
pub fn euler_slice<T, F>(mut f: F, t: T, y: &mut [T], h: T, work: &mut [T]) -> Result<(), LinalgError>
where
    T: Real,
    F: FnMut(T, &[T], &mut [T]),
{
    let Some(slope) = work.get_mut(..y.len()) else {
        return Err(LinalgError::DimensionMismatch);
    };
    f(t, y, slope);
    for (value, &k) in y.iter_mut().zip(slope.iter()) {
        *value += h * k;
    }
    Ok(())
}

/// Ein Schritt des klassischen Runge-Kutta-Verfahrens 4. Ordnung an Ort und
/// Stelle.
///
/// `f(t, y, dy)` schreibt die Ableitung nach `dy`. `work` braucht mindestens
/// `3 · y.len()` Elemente.
pub fn rk4_slice<T, F>(mut f: F, t: T, y: &mut [T], h: T, work: &mut [T]) -> Result<(), LinalgError>
where
    T: Real,
    F: FnMut(T, &[T], &mut [T]),
{
    let size = y.len();
    let Some(work) = size.checked_mul(3).and_then(|needed| work.get_mut(..needed)) else {
        return Err(LinalgError::DimensionMismatch);
    };
    let (slope, rest) = work.split_at_mut(size);
    let (stage, sum) = rest.split_at_mut(size);
    let half = h * T::from_f64(0.5);
    let two = T::from_f64(2.0);

    // Stufen nacheinander; `sum` sammelt k1 + 2·k2 + 2·k3 + k4.
    f(t, y, slope);
    sum.copy_from_slice(slope);
    for (weight, time, factor) in [(two, t + half, half), (two, t + half, half), (T::ONE, t + h, h)] {
        for ((target, &value), &k) in stage.iter_mut().zip(y.iter()).zip(slope.iter()) {
            *target = value + factor * k;
        }
        f(time, stage, slope);
        for (total, &k) in sum.iter_mut().zip(slope.iter()) {
            *total += weight * k;
        }
    }
    let scale = h / T::from_f64(6.0);
    for (value, &total) in y.iter_mut().zip(sum.iter()) {
        *value += scale * total;
    }
    Ok(())
}
