//! Dormand-Prince-Verfahren 5(4) mit Schrittweitensteuerung.
//!
//! Koeffizienten nach J. R. Dormand, P. J. Prince: „A family of embedded
//! Runge-Kutta formulae“, Journal of Computational and Applied Mathematics
//! 6 (1980), abgeglichen mit Wikipedia, „Dormand–Prince method“.

use super::error::{OdeError, OdeResult, OdeState};
use super::fixed::is_finite;
use crate::linalg::Vector;
use crate::scalar::RealElementary;

/// Einstellungen der Schrittweitensteuerung.
///
/// Ein Schritt wird angenommen, wenn der geschätzte Fehler jeder Komponente
/// im quadratischen Mittel höchstens `abs_tol + rel_tol·|y|` beträgt.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct OdeOptions<T> {
    /// Relative Toleranz.
    pub rel_tol: T,
    /// Absolute Toleranz.
    pub abs_tol: T,
    /// Erste Schrittweite; ohne Angabe ein Hundertstel des Intervalls.
    pub initial_step: Option<T>,
    /// Kleinste erlaubte Schrittweite.
    pub min_step: T,
    /// Höchstzahl der Schritte, angenommene und verworfene zusammen.
    pub max_steps: u32,
}

impl<T: RealElementary> Default for OdeOptions<T> {
    fn default() -> Self {
        Self {
            rel_tol: T::EPSILON.sqrt(),
            abs_tol: T::EPSILON,
            initial_step: None,
            min_step: T::ZERO,
            max_steps: 100_000,
        }
    }
}

impl<T: RealElementary> OdeOptions<T> {
    /// Setzt die relative Toleranz.
    #[inline]
    pub fn with_rel_tol(mut self, rel_tol: T) -> Self {
        self.rel_tol = rel_tol;
        self
    }

    /// Setzt die absolute Toleranz.
    #[inline]
    pub fn with_abs_tol(mut self, abs_tol: T) -> Self {
        self.abs_tol = abs_tol;
        self
    }

    /// Setzt die erste Schrittweite.
    #[inline]
    pub fn with_initial_step(mut self, initial_step: T) -> Self {
        self.initial_step = Some(initial_step);
        self
    }

    /// Setzt die kleinste erlaubte Schrittweite.
    #[inline]
    pub fn with_min_step(mut self, min_step: T) -> Self {
        self.min_step = min_step;
        self
    }

    /// Setzt die Höchstzahl der Schritte.
    #[inline]
    pub fn with_max_steps(mut self, max_steps: u32) -> Self {
        self.max_steps = max_steps;
        self
    }
}

/// Summe der Stufen, jede mit ihrem Koeffizienten gewichtet.
fn combine<T: RealElementary, const N: usize>(terms: &[(Vector<T, N>, T)]) -> Vector<T, N> {
    terms
        .iter()
        .fold(Vector::ZERO, |sum, &(stage, weight)| sum + stage * weight)
}

/// Löst `y' = f(t, y)` von `t0` nach `t1` mit dem Dormand-Prince-Verfahren.
///
/// Die Schrittweite passt sich an, sodass der geschätzte Fehler je Schritt
/// die Toleranzen aus `options` einhält. `t1 < t0` integriert rückwärts.
pub fn dormand_prince<T, F, const N: usize>(
    f: F,
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    options: &OdeOptions<T>,
) -> OdeResult<T, N>
where
    T: RealElementary,
    F: FnMut(T, Vector<T, N>) -> Vector<T, N>,
{
    dormand_prince_observed(f, t0, y0, t1, options, |_, _| {})
}

/// Wie [`dormand_prince`]; `observer` wird nach jedem angenommenen Schritt
/// mit Zeit und Lösung aufgerufen.
pub fn dormand_prince_observed<T, F, O, const N: usize>(
    mut f: F,
    t0: T,
    y0: Vector<T, N>,
    t1: T,
    options: &OdeOptions<T>,
    mut observer: O,
) -> OdeResult<T, N>
where
    T: RealElementary,
    F: FnMut(T, Vector<T, N>) -> Vector<T, N>,
    O: FnMut(T, Vector<T, N>),
{
    let tolerances_valid = options.rel_tol >= T::ZERO
        && options.abs_tol >= T::ZERO
        && (options.rel_tol > T::ZERO || options.abs_tol > T::ZERO);
    if !t0.is_finite() || !t1.is_finite() || !tolerances_valid || !is_finite(&y0) {
        return Err(OdeError::InvalidInput);
    }
    let mut state = OdeState { t: t0, y: y0, steps: 0, rejected: 0 };
    let span = t1 - t0;
    if span == T::ZERO {
        return Ok(state);
    }
    let direction = if span > T::ZERO { T::ONE } else { -T::ONE };
    let mut size = match options.initial_step {
        Some(step) => step.abs(),
        None => span.abs() * T::from_f64(0.01),
    };
    if size <= T::ZERO || size.is_nan() {
        return Err(OdeError::InvalidInput);
    }

    // Butcher-Tableau als Brüche; die Division geschieht zur Laufzeit.
    let r = |numerator: f64, denominator: f64| T::from_f64(numerator) / T::from_f64(denominator);
    let (c2, c3, c4, c5) = (r(1.0, 5.0), r(3.0, 10.0), r(4.0, 5.0), r(8.0, 9.0));
    let a21 = r(1.0, 5.0);
    let (a31, a32) = (r(3.0, 40.0), r(9.0, 40.0));
    let (a41, a42, a43) = (r(44.0, 45.0), r(-56.0, 15.0), r(32.0, 9.0));
    let (a51, a52, a53, a54) = (
        r(19372.0, 6561.0),
        r(-25360.0, 2187.0),
        r(64448.0, 6561.0),
        r(-212.0, 729.0),
    );
    let (a61, a62, a63, a64, a65) = (
        r(9017.0, 3168.0),
        r(-355.0, 33.0),
        r(46732.0, 5247.0),
        r(49.0, 176.0),
        r(-5103.0, 18656.0),
    );
    // Gewichte 5. Ordnung (zugleich letzte Zeile des Tableaus).
    let (b1, b3, b4, b5, b6) = (
        r(35.0, 384.0),
        r(500.0, 1113.0),
        r(125.0, 192.0),
        r(-2187.0, 6784.0),
        r(11.0, 84.0),
    );
    // Unterschied zu den Gewichten 4. Ordnung: liefert den Fehlerschätzer.
    let (e1, e3, e4, e5, e6, e7) = (
        b1 - r(5179.0, 57600.0),
        b3 - r(7571.0, 16695.0),
        b4 - r(393.0, 640.0),
        b5 - r(-92097.0, 339200.0),
        b6 - r(187.0, 2100.0),
        -r(1.0, 40.0),
    );
    let count = T::from_f64(f64::from(u32::try_from(N).unwrap_or(u32::MAX))).max(T::ONE);

    loop {
        if state.steps.saturating_add(state.rejected) >= options.max_steps {
            return Err(OdeError::TooManySteps(state));
        }
        let remaining = (t1 - state.t).abs();
        let last = size >= remaining;
        let taken = if last { remaining } else { size };
        let h = direction * taken;
        if taken < options.min_step || (!last && state.t + h == state.t) {
            return Err(OdeError::StepTooSmall(state));
        }

        let (t, y) = (state.t, state.y);
        let k1 = f(t, y);
        let k2 = f(t + c2 * h, y + combine(&[(k1, a21)]) * h);
        let k3 = f(t + c3 * h, y + combine(&[(k1, a31), (k2, a32)]) * h);
        let k4 = f(t + c4 * h, y + combine(&[(k1, a41), (k2, a42), (k3, a43)]) * h);
        let k5 = f(
            t + c5 * h,
            y + combine(&[(k1, a51), (k2, a52), (k3, a53), (k4, a54)]) * h,
        );
        let k6 = f(
            t + h,
            y + combine(&[(k1, a61), (k2, a62), (k3, a63), (k4, a64), (k5, a65)]) * h,
        );
        let next = y + combine(&[(k1, b1), (k3, b3), (k4, b4), (k5, b5), (k6, b6)]) * h;
        let k7 = f(t + h, next);
        let error =
            combine(&[(k1, e1), (k3, e3), (k4, e4), (k5, e5), (k6, e6), (k7, e7)]) * h;

        // Fehler im quadratischen Mittel, gemessen an der Toleranz.
        let mut sum = T::ZERO;
        for ((&old, &new), &deviation) in y.iter().zip(next.iter()).zip(error.iter()) {
            let scale = options.abs_tol + options.rel_tol * old.abs().max(new.abs());
            let ratio = deviation / scale;
            sum += ratio * ratio;
        }
        let measure = (sum / count).sqrt();
        let acceptable = is_finite(&next) && measure <= T::ONE;

        if acceptable {
            let reached = if last { t1 } else { t + h };
            state = OdeState { t: reached, y: next, steps: state.steps + 1, rejected: state.rejected };
            observer(state.t, state.y);
            if last {
                return Ok(state);
            }
        } else {
            state.rejected += 1;
        }

        // Neue Schrittweite aus dem Fehler: Faktor 0,9 · Fehler^(−1/5),
        // begrenzt auf [0,2; 5].
        let factor = if measure.is_nan() {
            T::from_f64(0.2)
        } else if measure == T::ZERO {
            T::from_f64(5.0)
        } else {
            (T::from_f64(0.9) * measure.powf(T::from_f64(-0.2)))
                .max(T::from_f64(0.2))
                .min(T::from_f64(5.0))
        };
        size = taken * factor;
    }
}
