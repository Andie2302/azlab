//! Newton-Verfahren für eine Unbekannte.

use super::error::{ScalarResult, SolveError};
use super::solution::Solution;
use super::tolerance::Tolerance;
use super::util::opposite_signs;
use crate::number::Dual;
use crate::scalar::Real;

/// Newton-Verfahren ab dem Startwert `x0`.
///
/// `f` liefert das Paar `(Funktionswert, Ableitung)`. Konvergiert sehr
/// schnell, aber nur in der Nähe einer Nullstelle; mit bekannter Klammer ist
/// [`newton_bracketed`] sicherer.
pub fn newton<T, F>(mut f: F, x0: T, tol: &Tolerance<T>) -> ScalarResult<T>
where
    T: Real,
    F: FnMut(T) -> (T, T),
{
    let mut x = x0;
    let (mut fx, mut slope) = f(x);
    for iteration in 1..=tol.max_iterations {
        let state = Solution { x, residual: fx.abs(), iterations: iteration - 1 };
        if fx == T::ZERO || tol.f_converged(fx) {
            return Ok(state);
        }
        if slope == T::ZERO {
            return Err(SolveError::Singular(state));
        }
        // Eine unendliche Ableitung ergäbe den Schritt null und sähe wie
        // Konvergenz aus.
        if !slope.is_finite() || !fx.is_finite() {
            return Err(SolveError::NonFinite(state));
        }
        let step = fx / slope;
        x -= step;
        (fx, slope) = f(x);
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

/// Newton-Verfahren mit automatischer Ableitung.
///
/// `f` wird mit Dualzahlen aufgerufen; die Ableitung entsteht dabei von
/// selbst. Man schreibt `f` einmal generisch über die Skalar-Traits.
pub fn newton_auto<T, F>(mut f: F, x0: T, tol: &Tolerance<T>) -> ScalarResult<T>
where
    T: Real,
    F: FnMut(Dual<T>) -> Dual<T>,
{
    newton(
        |x| {
            let value = f(Dual::variable(x));
            (value.re, value.eps)
        },
        x0,
        tol,
    )
}

/// Newton-Verfahren mit Klammer: Verlässt ein Newton-Schritt das Intervall
/// oder kommt er zu langsam voran, wird stattdessen halbiert.
///
/// `f` liefert das Paar `(Funktionswert, Ableitung)`. Konvergiert immer,
/// wenn `f(a)` und `f(b)` verschiedene Vorzeichen haben und `f` stetig ist.
pub fn newton_bracketed<T, F>(mut f: F, a: T, b: T, tol: &Tolerance<T>) -> ScalarResult<T>
where
    T: Real,
    F: FnMut(T) -> (T, T),
{
    if !(a.is_finite() && b.is_finite()) {
        return Err(SolveError::InvalidInput);
    }
    let (f_a, _) = f(a);
    let (f_b, _) = f(b);
    if f_a == T::ZERO {
        return Ok(Solution { x: a, residual: T::ZERO, iterations: 0 });
    }
    if f_b == T::ZERO {
        return Ok(Solution { x: b, residual: T::ZERO, iterations: 0 });
    }
    if !opposite_signs(f_a, f_b) {
        return Err(SolveError::NoBracket);
    }

    let half = T::from_f64(0.5);
    let two = T::from_f64(2.0);
    // negative: Seite mit f < 0, positive: Seite mit f > 0.
    let (mut negative, mut positive) = if f_a < T::ZERO { (a, b) } else { (b, a) };
    // Getrennt halbieren: Summen und Differenzen der Grenzen könnten
    // überlaufen.
    let mut x = half * a + half * b;
    let mut previous_step = (b - a).abs();
    let mut step = previous_step;
    let (mut fx, mut slope) = f(x);

    for iteration in 1..=tol.max_iterations {
        if fx == T::ZERO || tol.f_converged(fx) {
            return Ok(Solution { x, residual: fx.abs(), iterations: iteration - 1 });
        }
        // Liegt der Newton-Schritt echt innerhalb der Klammer? NaN und
        // unendliche Werte gelten als außerhalb.
        let candidate = x - fx / slope;
        let leaves_bracket =
            !(candidate > negative.min(positive) && candidate < negative.max(positive));
        let too_slow = (two * fx).abs() > (previous_step * slope).abs();
        previous_step = step;
        if leaves_bracket || too_slow || slope == T::ZERO || !slope.is_finite() {
            step = half * positive - half * negative;
            x = negative + step;
        } else {
            step = fx / slope;
            x -= step;
        }
        (fx, slope) = f(x);
        let state = Solution { x, residual: fx.abs(), iterations: iteration };
        if fx.is_nan() {
            return Err(SolveError::NonFinite(state));
        }
        if tol.x_converged(step, x) {
            return Ok(state);
        }
        if fx < T::ZERO {
            negative = x;
        } else {
            positive = x;
        }
    }
    Err(SolveError::NotConverged(Solution {
        x,
        residual: fx.abs(),
        iterations: tol.max_iterations,
    }))
}
