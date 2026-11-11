//! Skalar-Traits für [`Interval`].
//!
//! Die Rechengesetze der Traits gelten für Intervalle nur als Einschluss.

use super::repr::Interval;
use crate::scalar::{Algebra, Commutative, Field, Real, Ring};

impl<T: Real> Algebra for Interval<T> {
    const ZERO: Self = Self { lo: T::ZERO, hi: T::ZERO };
    const ONE: Self = Self { lo: T::ONE, hi: T::ONE };
}

impl<T: Real> Ring for Interval<T> {}

impl<T: Real> Commutative for Interval<T> {}

impl<T: Real> Field for Interval<T> {
    /// Ist `value` in `T` exakt darstellbar, entsteht ein Punktintervall,
    /// sonst ein Intervall, das `value` sicher enthält.
    fn from_f64(value: f64) -> Self {
        let rounded = T::from_f64(value);
        if rounded.to_f64() == value {
            Self::point(rounded)
        } else {
            Self::outward(rounded, rounded)
        }
    }
}
