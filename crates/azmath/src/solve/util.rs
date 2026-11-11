//! Kleine Hilfen für die Löser.

use crate::scalar::Real;

/// `true`, wenn `a` und `b` echt verschiedene Vorzeichen haben.
/// Null und NaN zählen nicht als Vorzeichen.
#[inline]
pub(super) fn opposite_signs<T: Real>(a: T, b: T) -> bool {
    (a < T::ZERO && b > T::ZERO) || (a > T::ZERO && b < T::ZERO)
}
