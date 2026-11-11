//! Nelder-Mead-Verfahren (Simplex-Verfahren ohne Ableitungen).

use crate::linalg::Vector;
use crate::scalar::Real;
use crate::solve::{Solution, SolveError, SystemResult, Tolerance};

/// Startsimplex um `x0`: `x0` selbst und je ein Punkt, der in einer
/// Koordinatenrichtung um `step` verschoben ist.
///
/// `P` muss `N + 1` sein, sonst `None`.
pub fn simplex_around<T: Real, const N: usize, const P: usize>(
    x0: Vector<T, N>,
    step: T,
) -> Option<[Vector<T, N>; P]> {
    if P != N + 1 {
        return None;
    }
    Some(core::array::from_fn(|point| {
        let mut corner = x0;
        if let Some(component) = point.checked_sub(1).and_then(|axis| corner.get_mut(axis)) {
            *component += step;
        }
        corner
    }))
}

/// Minimiert `f` mit dem Nelder-Mead-Verfahren, ausgehend vom Simplex
/// `simplex` mit `P = N + 1` Ecken.
///
/// Braucht nur Funktionswerte, keine Ableitungen, und kommt deshalb auch mit
/// nicht glatten Funktionen zurecht. Konvergiert langsamer als die
/// gradientenbasierten Verfahren und nur für wenige Unbekannte gut.
///
/// Abbruch, wenn die Spanne der Funktionswerte an den Ecken höchstens
/// `abs_f + rel_x · |bester Wert|` ist. Im Ergebnis ist `x` die beste Ecke
/// und `residual` diese Spanne.
pub fn nelder_mead<T, F, const N: usize, const P: usize>(
    mut f: F,
    simplex: [Vector<T, N>; P],
    tol: &Tolerance<T>,
) -> SystemResult<T, N>
where
    T: Real,
    F: FnMut(Vector<T, N>) -> T,
{
    if P != N + 1 || N == 0 {
        return Err(SolveError::InvalidInput);
    }
    let half = T::from_f64(0.5);
    let two = T::from_f64(2.0);
    let mut corners = simplex.map(|x| (f(x), x));

    for iteration in 0..=tol.max_iterations {
        // Aufsteigend nach Funktionswert sortieren (NaN ans Ende).
        for sorted in 1..P {
            let mut slot = sorted;
            while slot > 0
                && (corners[slot].0 < corners[slot - 1].0 || corners[slot - 1].0.is_nan())
            {
                corners.swap(slot - 1, slot);
                slot -= 1;
            }
        }
        let (best, worst) = (corners[0], corners[P - 1]);
        let second_worst = corners[P - 2].0;
        let spread = (worst.0 - best.0).abs();
        let state = Solution { x: best.1, residual: spread, iterations: iteration };
        if !best.0.is_finite() {
            return Err(SolveError::NonFinite(state));
        }
        if spread <= tol.abs_f + tol.rel_x * best.0.abs() {
            return Ok(state);
        }
        if iteration == tol.max_iterations {
            return Err(SolveError::NotConverged(state));
        }

        // Schwerpunkt aller Ecken außer der schlechtesten.
        let sum = corners[..P - 1].iter().fold(Vector::ZERO, |acc, corner| acc + corner.1);
        let centroid = sum * T::from_f64(N as f64).recip();
        let along = |factor: T| centroid + (centroid - worst.1) * factor;

        let reflected = along(T::ONE);
        let f_reflected = f(reflected);
        if f_reflected < best.0 {
            // Sehr gut: noch weiter in diese Richtung versuchen.
            let expanded = along(two);
            let f_expanded = f(expanded);
            corners[P - 1] = if f_expanded < f_reflected {
                (f_expanded, expanded)
            } else {
                (f_reflected, reflected)
            };
        } else if f_reflected < second_worst {
            corners[P - 1] = (f_reflected, reflected);
        } else {
            // Schlecht: zusammenziehen, außen oder innen.
            let contracted = if f_reflected < worst.0 { along(half) } else { along(-half) };
            let f_contracted = f(contracted);
            if f_contracted < worst.0.min(f_reflected) {
                corners[P - 1] = (f_contracted, contracted);
            } else {
                // Auch das half nicht: alle Ecken zur besten hin schrumpfen.
                for corner in corners.iter_mut().skip(1) {
                    let x = best.1 + (corner.1 - best.1) * half;
                    *corner = (f(x), x);
                }
            }
        }
    }
    Err(SolveError::InvalidInput)
}
