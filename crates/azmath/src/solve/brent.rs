//! Verfahren von Brent.

use super::error::{ScalarResult, SolveError};
use super::solution::Solution;
use super::tolerance::Tolerance;
use super::util::opposite_signs;
use crate::scalar::Real;

/// Sucht eine Nullstelle von `f` zwischen `a` und `b` mit dem Verfahren von
/// Brent.
///
/// Verbindet die Sicherheit der Bisektion mit der Geschwindigkeit von
/// Sekante und inverser quadratischer Interpolation. Die Näherung verlässt
/// die Klammer nie. Das Standardverfahren, wenn eine Klammer bekannt ist.
pub fn brent<T, F>(mut f: F, a: T, b: T, tol: &Tolerance<T>) -> ScalarResult<T>
where
    T: Real,
    F: FnMut(T) -> T,
{
    if !(a.is_finite() && b.is_finite()) {
        return Err(SolveError::InvalidInput);
    }
    let (mut a, mut b) = (a, b);
    let mut fa = f(a);
    let mut fb = f(b);
    if fa == T::ZERO {
        return Ok(Solution { x: a, residual: T::ZERO, iterations: 0 });
    }
    if fb == T::ZERO {
        return Ok(Solution { x: b, residual: T::ZERO, iterations: 0 });
    }
    if !opposite_signs(fa, fb) {
        return Err(SolveError::NoBracket);
    }

    let half = T::from_f64(0.5);
    let two = T::from_f64(2.0);
    let three = T::from_f64(3.0);
    // b ist die beste Näherung, c der Gegenpunkt der Klammer, a der
    // vorherige Wert von b.
    let mut c = b;
    let mut fc = fb;
    let mut d = b - a;
    let mut e = d;

    for iteration in 1..=tol.max_iterations {
        if (fb > T::ZERO && fc > T::ZERO) || (fb < T::ZERO && fc < T::ZERO) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        let tol1 = two * T::EPSILON * b.abs() + half * (tol.abs_x + tol.rel_x * b.abs());
        // Getrennt halbieren: `c - b` könnte überlaufen.
        let xm = half * c - half * b;
        if xm.abs() <= tol1 || fb == T::ZERO || tol.f_converged(fb) {
            let state = Solution { x: b, residual: fb.abs(), iterations: iteration - 1 };
            // Kein Erfolg mit unendlichem Restfehler (z. B. an einer Polstelle).
            return if fb.is_finite() { Ok(state) } else { Err(SolveError::NonFinite(state)) };
        }

        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            // Interpolation versuchen.
            let s = fb / fa;
            let (mut p, mut q) = if a == c {
                (two * xm * s, T::ONE - s)
            } else {
                let r0 = fa / fc;
                let r = fb / fc;
                (
                    s * (two * xm * r0 * (r0 - r) - (b - a) * (r - T::ONE)),
                    (r0 - T::ONE) * (r - T::ONE) * (s - T::ONE),
                )
            };
            if p > T::ZERO {
                q = -q;
            }
            p = p.abs();
            let bound1 = three * xm * q - (tol1 * q).abs();
            let bound2 = (e * q).abs();
            if two * p < bound1.min(bound2) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            // Zu langsam: Bisektion.
            d = xm;
            e = d;
        }

        a = b;
        fa = fb;
        if d.abs() > tol1 {
            b += d;
        } else {
            b += tol1.copysign(xm);
        }
        fb = f(b);
        if fb.is_nan() {
            return Err(SolveError::NonFinite(Solution {
                x: b,
                residual: fb,
                iterations: iteration,
            }));
        }
    }
    Err(SolveError::NotConverged(Solution {
        x: b,
        residual: fb.abs(),
        iterations: tol.max_iterations,
    }))
}
