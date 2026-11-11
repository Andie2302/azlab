//! Bisektion.

use super::error::{ScalarResult, SolveError};
use super::solution::Solution;
use super::tolerance::Tolerance;
use super::util::opposite_signs;
use crate::scalar::Real;

/// Sucht eine Nullstelle von `f` zwischen `a` und `b` durch fortgesetztes
/// Halbieren.
///
/// Langsam, aber unbedingt zuverlässig: Haben `f(a)` und `f(b)` verschiedene
/// Vorzeichen und ist `f` stetig, wird eine Nullstelle gefunden. Für die
/// meisten Zwecke ist [`brent`](super::brent()) die bessere Wahl.
pub fn bisection<T, F>(mut f: F, a: T, b: T, tol: &Tolerance<T>) -> ScalarResult<T>
where
    T: Real,
    F: FnMut(T) -> T,
{
    if !(a.is_finite() && b.is_finite()) {
        return Err(SolveError::InvalidInput);
    }
    let (mut lo, mut hi) = if a <= b { (a, b) } else { (b, a) };
    let mut f_lo = f(lo);
    let f_hi = f(hi);
    if f_lo == T::ZERO {
        return Ok(Solution { x: lo, residual: T::ZERO, iterations: 0 });
    }
    if f_hi == T::ZERO {
        return Ok(Solution { x: hi, residual: T::ZERO, iterations: 0 });
    }
    if !opposite_signs(f_lo, f_hi) {
        return Err(SolveError::NoBracket);
    }

    let half = T::from_f64(0.5);
    let mut x = lo;
    let mut fx = f_lo;
    for iteration in 1..=tol.max_iterations {
        // Getrennt halbieren: `hi - lo` könnte überlaufen.
        x = lo * half + hi * half;
        fx = f(x);
        let state = Solution { x, residual: fx.abs(), iterations: iteration };
        if fx.is_nan() {
            return Err(SolveError::NonFinite(state));
        }
        // `x == lo || x == hi`: das Intervall lässt sich nicht weiter teilen.
        if fx == T::ZERO
            || tol.f_converged(fx)
            || x == lo
            || x == hi
            || tol.x_converged(hi - lo, x)
        {
            return Ok(state);
        }
        if opposite_signs(f_lo, fx) {
            hi = x;
        } else {
            lo = x;
            f_lo = fx;
        }
    }
    Err(SolveError::NotConverged(Solution {
        x,
        residual: fx.abs(),
        iterations: tol.max_iterations,
    }))
}
