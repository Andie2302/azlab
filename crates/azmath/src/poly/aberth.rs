//! Alle komplexen Nullstellen eines Polynoms mit dem Aberth-Ehrlich-Verfahren.

use super::eval::horner_with_derivative;
use crate::number::Complex;
use crate::scalar::RealElementary;
use crate::solve::{Solution, SolveError, Tolerance};

/// Berechnet alle komplexen Nullstellen von
/// `c[0] + c[1]·x + … + c[n]·xⁿ` gleichzeitig.
///
/// `roots` muss genau so lang sein wie der Grad (höchster Koeffizient
/// ungleich null). Rückgabe: Zahl der Iterationen.
///
/// Das Verfahren ist schnell und funktioniert für jeden Grad, ist aber ein
/// Näherungsverfahren **ohne Garantie**. Wo keine reelle Nullstelle fehlen
/// darf, ist [`isolate_real_roots`](super::isolate_real_roots) das richtige
/// Werkzeug.
///
/// Fehler: [`SolveError::InvalidInput`], wenn alle Koeffizienten null sind
/// oder `roots` zu lang ist; [`SolveError::CapacityExceeded`], wenn `roots`
/// zu kurz ist. Bei `NotConverged` stehen die letzten Näherungen in `roots`.
// Die Schleife über den Index braucht gleichzeitig Lese- und Schreibzugriff.
#[allow(clippy::needless_range_loop)]
pub fn complex_roots<T: RealElementary>(
    coefficients: &[Complex<T>],
    roots: &mut [Complex<T>],
    tol: &Tolerance<T>,
) -> Result<u32, SolveError<(), T>> {
    let zero = Complex::new(T::ZERO, T::ZERO);
    let one = Complex::new(T::ONE, T::ZERO);

    let Some(degree) = coefficients.iter().rposition(|c| *c != zero) else {
        return Err(SolveError::InvalidInput);
    };
    let Some(coefficients) = coefficients.get(..=degree) else {
        return Err(SolveError::InvalidInput);
    };
    if roots.len() < degree {
        return Err(SolveError::CapacityExceeded);
    }
    if roots.len() > degree {
        return Err(SolveError::InvalidInput);
    }
    let Some((&lead, lower)) = coefficients.split_last() else {
        return Err(SolveError::InvalidInput);
    };
    if degree == 0 {
        return Ok(0);
    }

    // Startwerte: gleichmäßig auf einem Kreis, der alle Nullstellen enthält
    // (Schranke von Cauchy), leicht verdreht gegen Symmetrien.
    let mut radius = T::ZERO;
    let mut count = T::ZERO;
    for &c in lower {
        radius = radius.max((c / lead).norm());
        count += T::ONE;
    }
    radius += T::ONE;
    let step = T::from_f64(2.0) * T::PI / count;
    let mut angle = T::from_f64(0.4);
    for root in roots.iter_mut() {
        *root = Complex::from_polar(radius, angle);
        angle += step;
    }

    let mut residual = T::ZERO;
    for iteration in 1..=tol.max_iterations {
        let mut converged = true;
        residual = T::ZERO;
        for i in 0..degree {
            let z = roots[i];
            let (value, slope) = horner_with_derivative(coefficients, z);
            let magnitude = value.norm();
            residual = residual.max(magnitude);
            if magnitude == T::ZERO {
                continue;
            }
            // Newton-Schritt, korrigiert um die Abstoßung der anderen
            // Näherungen.
            let newton = value / slope;
            let mut repulsion = zero;
            for (j, &other) in roots.iter().enumerate() {
                if j != i {
                    repulsion += one / (z - other);
                }
            }
            let correction = newton / (one - newton * repulsion);
            if correction.re.is_nan() || correction.im.is_nan() {
                return Err(SolveError::NonFinite(Solution {
                    x: (),
                    residual,
                    iterations: iteration,
                }));
            }
            roots[i] = z - correction;
            if !tol.x_converged(correction.norm(), z.norm()) {
                converged = false;
            }
        }
        if converged {
            return Ok(iteration);
        }
    }
    Err(SolveError::NotConverged(Solution {
        x: (),
        residual,
        iterations: tol.max_iterations,
    }))
}
