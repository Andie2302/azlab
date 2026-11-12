//! Kleine Hilfen für die Zerlegungen.

use crate::scalar::{Conjugate, Field, Real};

/// `true` nur für echt positive Werte; NaN gilt als nicht positiv.
#[inline]
pub(super) fn is_positive<R: Real>(value: R) -> bool {
    value > R::ZERO
}

/// Berechnet `divisor⁻¹ · value`, ohne den Kehrwert von `divisor` selbst zu
/// bilden.
///
/// Der Kehrwert eines sehr kleinen Diagonalelements läuft über, auch wenn
/// das Ergebnis darstellbar ist (etwa `1e-320 / 1e-320`). Deshalb werden
/// beide Seiten zuerst durch den Betrag des Teilers geteilt. Die Reihenfolge
/// bleibt erhalten, das Ergebnis stimmt also auch für nicht kommutative
/// Typen wie Quaternionen.
#[inline]
pub(crate) fn left_divide<T: Field + Conjugate>(divisor: T, value: T) -> T {
    let scale = divisor.max_abs();
    if is_positive(scale) && scale.is_finite() {
        divisor.div_real(scale).recip() * value.div_real(scale)
    } else {
        divisor.recip() * value
    }
}
