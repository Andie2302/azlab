//! Lineare Interpolation.

use crate::scalar::Real;

/// Lineare Interpolation zwischen `a` (bei `t = 0`) und `b` (bei `t = 1`).
#[inline]
pub fn lerp<T: Real>(a: T, b: T, t: T) -> T {
    a + (b - a) * t
}

/// Index `i` des Teilintervalls `[xs[i], xs[i + 1]]`, das `x` enthält.
///
/// `None`, wenn es weniger als zwei Stützstellen gibt oder `x` außerhalb
/// liegt (auch für NaN).
pub(super) fn locate<T: Real>(xs: &[T], x: T) -> Option<usize> {
    let (first, last) = (*xs.first()?, *xs.last()?);
    if xs.len() < 2 || !(first <= x && x <= last) {
        return None;
    }
    let upper = xs.partition_point(|&value| value <= x);
    Some(upper.saturating_sub(1).min(xs.len() - 2))
}

/// Stückweise lineare Interpolation der Punkte `(xs[i], ys[i])` an der
/// Stelle `x`.
///
/// `xs` muss aufsteigend sortiert sein. `None`, wenn die Längen verschieden
/// sind, es weniger als zwei Punkte gibt oder `x` außerhalb liegt.
pub fn linear<T: Real>(xs: &[T], ys: &[T], x: T) -> Option<T> {
    if xs.len() != ys.len() {
        return None;
    }
    let index = locate(xs, x)?;
    let (x0, x1) = (*xs.get(index)?, *xs.get(index + 1)?);
    let (y0, y1) = (*ys.get(index)?, *ys.get(index + 1)?);
    if x1 == x0 {
        return Some(y0);
    }
    Some(lerp(y0, y1, (x - x0) / (x1 - x0)))
}
