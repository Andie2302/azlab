//! Explizite Verfahren mit fester Schrittweite.

use super::error::{OdeError, OdeResult, OdeState};
use crate::linalg::Vector;
use crate::scalar::Real;

/// `true`, wenn alle Komponenten endlich sind.
pub(super) fn is_finite<T: Real, const N: usize>(y: &Vector<T, N>) -> bool {
    y.iter().all(|component| component.is_finite())
}

/// Ein Schritt des expliziten Euler-Verfahrens (1. Ordnung).
#[inline]
pub fn euler_step<T, F, const N: usize>(f: &mut F, t: T, y: Vector<T, N>, h: T) -> Vector<T, N>
where
    T: Real,
    F: FnMut(T, Vector<T, N>) -> Vector<T, N>,
{
    y + f(t, y) * h
}

/// Ein Schritt des klassischen Runge-Kutta-Verfahrens (4. Ordnung).
pub fn rk4_step<T, F, const N: usize>(f: &mut F, t: T, y: Vector<T, N>, h: T) -> Vector<T, N>
where
    T: Real,
    F: FnMut(T, Vector<T, N>) -> Vector<T, N>,
{
    let half = h * T::from_f64(0.5);
    let two = T::from_f64(2.0);
    let k1 = f(t, y);
    let k2 = f(t + half, y + k1 * half);
    let k3 = f(t + half, y + k2 * half);
    let k4 = f(t + h, y + k3 * h);
    y + (k1 + k2 * two + k3 * two + k4) * (h / T::from_f64(6.0))
}

/// Führt `steps` gleich große Schritte von `t0` nach `t1` aus.
pub(super) fn fixed_steps<T, S, const N: usize>(
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    steps: u32,
    mut step: S,
) -> OdeResult<T, N>
where
    T: Real,
    S: FnMut(T, Vector<T, N>, T, OdeState<T, N>) -> Result<Vector<T, N>, OdeError<T, N>>,
{
    if steps == 0 || !t0.is_finite() || !t1.is_finite() {
        return Err(OdeError::InvalidInput);
    }
    let h = (t1 - t0) / T::from_f64(f64::from(steps));
    let mut state = OdeState { t: t0, y: y0, steps: 0, rejected: 0 };
    for index in 1..=steps {
        let next = step(state.t, state.y, h, state)?;
        if !is_finite(&next) {
            return Err(OdeError::NonFinite(state));
        }
        // Die Zeit wird aus dem Zähler berechnet, damit sich Rundungsfehler
        // nicht aufsummieren; der letzte Schritt trifft t1 exakt.
        let t = if index == steps { t1 } else { t0 + h * T::from_f64(f64::from(index)) };
        state = OdeState { t, y: next, steps: index, rejected: 0 };
    }
    Ok(state)
}

/// Explizites Euler-Verfahren von `t0` nach `t1` in `steps` Schritten.
///
/// Einfach, aber nur 1. Ordnung; für genaue Ergebnisse [`rk4`] oder
/// [`dormand_prince`](super::dormand_prince) nehmen.
pub fn euler<T, F, const N: usize>(
    mut f: F,
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    steps: u32,
) -> OdeResult<T, N>
where
    T: Real,
    F: FnMut(T, Vector<T, N>) -> Vector<T, N>,
{
    fixed_steps(t0, y0, t1, steps, |t, y, h, _| Ok(euler_step(&mut f, t, y, h)))
}

/// Klassisches Runge-Kutta-Verfahren 4. Ordnung von `t0` nach `t1` in
/// `steps` Schritten.
pub fn rk4<T, F, const N: usize>(
    mut f: F,
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    steps: u32,
) -> OdeResult<T, N>
where
    T: Real,
    F: FnMut(T, Vector<T, N>) -> Vector<T, N>,
{
    fixed_steps(t0, y0, t1, steps, |t, y, h, _| Ok(rk4_step(&mut f, t, y, h)))
}
