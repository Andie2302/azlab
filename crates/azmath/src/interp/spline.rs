//! Natürliche kubische Splines.

use super::error::InterpError;
use super::linear::locate;
use crate::scalar::Real;

/// Berechnet die zweiten Ableitungen `m` des natürlichen kubischen Splines
/// durch die Punkte `(xs[i], ys[i])`.
///
/// „Natürlich“ heißt: Die zweite Ableitung ist an beiden Enden null. `work`
/// ist Arbeitsspeicher derselben Länge. Danach wertet [`spline_eval`] den
/// Spline aus.
///
/// Fehler, wenn die Längen verschieden sind, es weniger als zwei Punkte
/// gibt oder `xs` nicht streng aufsteigend ist.
// Die Rekursion braucht Nachbarn links und rechts in mehreren Feldern.
#[allow(clippy::needless_range_loop)]
pub fn natural_spline<T: Real>(
    xs: &[T],
    ys: &[T],
    m: &mut [T],
    work: &mut [T],
) -> Result<(), InterpError> {
    let count = xs.len();
    if ys.len() != count || m.len() != count || work.len() != count {
        return Err(InterpError::LengthMismatch);
    }
    if count < 2 {
        return Err(InterpError::TooFewPoints);
    }
    if !xs.windows(2).all(|pair| pair[0] < pair[1]) {
        return Err(InterpError::NotIncreasing);
    }

    // Tridiagonales System, gelöst mit dem Thomas-Algorithmus.
    let two = T::from_f64(2.0);
    m[0] = T::ZERO;
    work[0] = T::ZERO;
    for i in 1..count - 1 {
        let ratio = (xs[i] - xs[i - 1]) / (xs[i + 1] - xs[i - 1]);
        let pivot = ratio * m[i - 1] + two;
        m[i] = (ratio - T::ONE) / pivot;
        let slope_change = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
            - (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
        work[i] = (T::from_f64(6.0) * slope_change / (xs[i + 1] - xs[i - 1])
            - ratio * work[i - 1])
            / pivot;
    }
    m[count - 1] = T::ZERO;
    for i in (0..count - 1).rev() {
        m[i] = m[i] * m[i + 1] + work[i];
    }
    Ok(())
}

/// Stützwerte des Teilintervalls, das `x` enthält:
/// `(Breite, a, b, y0, y1, m0, m1)` mit den Gewichten `a`, `b` der Enden.
fn segment<T: Real>(xs: &[T], ys: &[T], m: &[T], x: T) -> Option<(T, T, T, T, T, T, T)> {
    if ys.len() != xs.len() || m.len() != xs.len() {
        return None;
    }
    let index = locate(xs, x)?;
    let (x0, x1) = (*xs.get(index)?, *xs.get(index + 1)?);
    let width = x1 - x0;
    Some((
        width,
        (x1 - x) / width,
        (x - x0) / width,
        *ys.get(index)?,
        *ys.get(index + 1)?,
        *m.get(index)?,
        *m.get(index + 1)?,
    ))
}

/// Wert des Splines an der Stelle `x`; `m` stammt aus [`natural_spline`].
///
/// `None`, wenn `x` außerhalb der Stützstellen liegt oder die Längen nicht
/// passen.
pub fn spline_eval<T: Real>(xs: &[T], ys: &[T], m: &[T], x: T) -> Option<T> {
    let (width, a, b, y0, y1, m0, m1) = segment(xs, ys, m, x)?;
    let curvature = ((a * a * a - a) * m0 + (b * b * b - b) * m1) * width * width
        / T::from_f64(6.0);
    Some(a * y0 + b * y1 + curvature)
}

/// Erste Ableitung des Splines an der Stelle `x`.
///
/// `None`, wenn `x` außerhalb der Stützstellen liegt oder die Längen nicht
/// passen.
pub fn spline_derivative<T: Real>(xs: &[T], ys: &[T], m: &[T], x: T) -> Option<T> {
    let (width, a, b, y0, y1, m0, m1) = segment(xs, ys, m, x)?;
    let three = T::from_f64(3.0);
    let bend = ((three * b * b - T::ONE) * m1 - (three * a * a - T::ONE) * m0) * width
        / T::from_f64(6.0);
    Some((y1 - y0) / width + bend)
}

/// Natürlicher kubischer Spline durch `N` Punkte, ohne Heap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubicSpline<T, const N: usize> {
    xs: [T; N],
    ys: [T; N],
    m: [T; N],
}

impl<T: Real, const N: usize> CubicSpline<T, N> {
    /// Baut den Spline durch die Punkte `(xs[i], ys[i])`.
    ///
    /// Fehler, wenn es weniger als zwei Punkte gibt oder `xs` nicht streng
    /// aufsteigend ist.
    pub fn natural(xs: [T; N], ys: [T; N]) -> Result<Self, InterpError> {
        let mut m = [T::ZERO; N];
        let mut work = [T::ZERO; N];
        natural_spline(&xs, &ys, &mut m, &mut work)?;
        Ok(Self { xs, ys, m })
    }

    /// Wert an der Stelle `x`; `None` außerhalb der Stützstellen.
    #[inline]
    pub fn eval(&self, x: T) -> Option<T> {
        spline_eval(&self.xs, &self.ys, &self.m, x)
    }

    /// Erste Ableitung an der Stelle `x`; `None` außerhalb der Stützstellen.
    #[inline]
    pub fn derivative(&self, x: T) -> Option<T> {
        spline_derivative(&self.xs, &self.ys, &self.m, x)
    }
}
