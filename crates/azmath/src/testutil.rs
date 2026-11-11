//! Hilfen für Tests.

use crate::scalar::Real;

/// `true`, wenn `a` und `b` bis auf `ulps` Maschinengenauigkeiten
/// übereinstimmen (relativ zu `max(1, |b|)`).
pub(crate) fn close<T: Real>(a: T, b: T, ulps: f64) -> bool {
    (a - b).abs() <= T::from_f64(ulps) * T::EPSILON * b.abs().max(T::ONE)
}
