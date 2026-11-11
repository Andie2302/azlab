//! Kleine Hilfen für die Zerlegungen.

use crate::scalar::Real;

/// `true` nur für echt positive Werte; NaN gilt als nicht positiv.
#[inline]
pub(super) fn is_positive<R: Real>(value: R) -> bool {
    value > R::ZERO
}
