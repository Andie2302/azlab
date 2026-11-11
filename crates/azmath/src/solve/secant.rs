//! Sekantenverfahren.

use super::error::{ScalarResult, SolveError};
use super::solution::Solution;
use super::tolerance::Tolerance;
use crate::scalar::Real;

/// Sucht eine Nullstelle von `f` mit dem Sekantenverfahren, ausgehend von
/// den beiden Startwerten `x0` und `x1`.
///
/// Braucht weder Ableitung noch Klammer, konvergiert aber nur in der Nähe
/// einer Nullstelle zuverlässig.
pub fn secant<T, F>(mut f: F, x0: T, x1: T, tol: &Tolerance<T>) -> ScalarResult<T>
where
    T: Real,
    F: FnMut(T) -> T,
{
    let mut previous = x0;
    let mut x = x1;
    let mut f_previous = f(previous);
    let mut fx = f(x);
    for iteration in 1..=tol.max_iterations {
        let state = Solution { x, residual: fx.abs(), iterations: iteration - 1 };
        if fx == T::ZERO || tol.f_converged(fx) {
            return Ok(state);
        }
        let rise = fx - f_previous;
        if rise == T::ZERO {
            return Err(SolveError::Singular(state));
        }
        let step = fx * (x - previous) / rise;
        previous = x;
        f_previous = fx;
        x -= step;
        fx = f(x);
        let state = Solution { x, residual: fx.abs(), iterations: iteration };
        if !x.is_finite() || fx.is_nan() {
            return Err(SolveError::NonFinite(state));
        }
        if tol.x_converged(step, x) {
            return Ok(state);
        }
    }
    Err(SolveError::NotConverged(Solution {
        x,
        residual: fx.abs(),
        iterations: tol.max_iterations,
    }))
}
