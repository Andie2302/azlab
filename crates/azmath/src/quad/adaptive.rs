//! Adaptive Simpson-Quadratur.

use crate::scalar::Real;
use core::fmt;

/// Tiefe des Stapels noch nicht bearbeiteter Teilintervalle.
const STACK: usize = 64;

/// Einstellungen der adaptiven Quadratur.
///
/// Angestrebt wird ein Gesamtfehler von höchstens
/// `abs_tol + rel_tol · |Integral|`.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct QuadOptions<T> {
    /// Absolute Toleranz.
    pub abs_tol: T,
    /// Relative Toleranz.
    pub rel_tol: T,
    /// Höchstzahl der Funktionsauswertungen.
    pub max_evaluations: u32,
}

impl<T: Real> Default for QuadOptions<T> {
    fn default() -> Self {
        Self {
            abs_tol: T::ZERO,
            rel_tol: T::EPSILON * T::from_f64(1024.0),
            max_evaluations: 100_000,
        }
    }
}

impl<T: Real> QuadOptions<T> {
    /// Setzt die absolute Toleranz.
    #[inline]
    pub fn with_abs_tol(mut self, abs_tol: T) -> Self {
        self.abs_tol = abs_tol;
        self
    }

    /// Setzt die relative Toleranz.
    #[inline]
    pub fn with_rel_tol(mut self, rel_tol: T) -> Self {
        self.rel_tol = rel_tol;
        self
    }

    /// Setzt die Höchstzahl der Funktionsauswertungen.
    #[inline]
    pub fn with_max_evaluations(mut self, max_evaluations: u32) -> Self {
        self.max_evaluations = max_evaluations;
        self
    }
}

/// Ergebnis einer Quadratur mit Fehlerschätzung.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Quadrature<T> {
    /// Näherung des Integrals.
    pub value: T,
    /// Geschätzter Betrag des Fehlers.
    pub error: T,
    /// Anzahl der Funktionsauswertungen.
    pub evaluations: u32,
}

/// Fehler der adaptiven Quadratur.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum QuadError<T> {
    /// Ungültige Eingabe, z. B. nicht endliche Grenzen.
    InvalidInput,
    /// Die Toleranz wurde nicht erreicht; enthält die beste Näherung.
    NotConverged(Quadrature<T>),
}

impl<T> fmt::Display for QuadError<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => "ungültige Eingabe",
            Self::NotConverged(_) => "Toleranz nicht erreicht",
        })
    }
}

impl<T: fmt::Debug> core::error::Error for QuadError<T> {}

/// Teilintervall mit bereits bekannten Funktionswerten.
#[derive(Clone, Copy)]
struct Panel<T> {
    a: T,
    fa: T,
    m: T,
    fm: T,
    b: T,
    fb: T,
    whole: T,
    tolerance: T,
}

/// Integral von `f` über [a, b] mit adaptiver Simpsonregel.
///
/// Das Intervall wird dort weiter geteilt, wo die Fehlerschätzung zu groß
/// ist. Arbeitet ohne Rekursion und ohne Heap.
pub fn adaptive_simpson<T, F>(
    mut f: F,
    a: T,
    b: T,
    options: &QuadOptions<T>,
) -> Result<Quadrature<T>, QuadError<T>>
where
    T: Real,
    F: FnMut(T) -> T,
{
    let tolerances_valid = options.abs_tol >= T::ZERO
        && options.rel_tol >= T::ZERO
        && (options.abs_tol > T::ZERO || options.rel_tol > T::ZERO);
    if !a.is_finite() || !b.is_finite() || !tolerances_valid {
        return Err(QuadError::InvalidInput);
    }
    let half = T::from_f64(0.5);
    let sixth = T::ONE / T::from_f64(6.0);
    let fifteen = T::from_f64(15.0);
    let rule = |a: T, fa: T, fm: T, b: T, fb: T| (b - a) * sixth * (fa + T::from_f64(4.0) * fm + fb);

    let m = (a + b) * half;
    let (fa, fm, fb) = (f(a), f(m), f(b));
    let whole = rule(a, fa, fm, b, fb);
    let mut evaluations = 3_u32;
    let tolerance = options.abs_tol + options.rel_tol * whole.abs();

    let first = Panel { a, fa, m, fm, b, fb, whole, tolerance };
    let mut stack = [first; STACK];
    let mut depth = 1;
    let mut value = T::ZERO;
    let mut error = T::ZERO;
    let mut converged = true;

    while depth > 0 {
        depth -= 1;
        let panel = stack[depth];
        let left_mid = (panel.a + panel.m) * half;
        let right_mid = (panel.m + panel.b) * half;
        let (f_left, f_right) = (f(left_mid), f(right_mid));
        evaluations = evaluations.saturating_add(2);
        let left = rule(panel.a, panel.fa, f_left, panel.m, panel.fm);
        let right = rule(panel.m, panel.fm, f_right, panel.b, panel.fb);
        let delta = left + right - panel.whole;

        let accurate = delta.abs() <= fifteen * panel.tolerance;
        // Nicht weiter teilbar: Mittelpunkte fallen mit den Rändern zusammen.
        let exhausted = left_mid <= panel.a || right_mid >= panel.b;
        let out_of_budget = depth + 2 > STACK || evaluations >= options.max_evaluations;
        if accurate || exhausted || out_of_budget {
            if !accurate && !exhausted {
                converged = false;
            }
            value += left + right + delta / fifteen;
            error += (delta / fifteen).abs();
        } else {
            let tolerance = panel.tolerance * half;
            stack[depth] = Panel {
                a: panel.m,
                fa: panel.fm,
                m: right_mid,
                fm: f_right,
                b: panel.b,
                fb: panel.fb,
                whole: right,
                tolerance,
            };
            stack[depth + 1] = Panel {
                a: panel.a,
                fa: panel.fa,
                m: left_mid,
                fm: f_left,
                b: panel.m,
                fb: panel.fm,
                whole: left,
                tolerance,
            };
            depth += 2;
        }
    }

    let result = Quadrature { value, error, evaluations };
    if converged && value.is_finite() {
        Ok(result)
    } else {
        Err(QuadError::NotConverged(result))
    }
}
