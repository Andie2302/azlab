//! Skalar-Traits für [`Quaternion`].
//!
//! Quaternionen bilden einen Schiefkörper: `Ring` und `Field`, aber nicht
//! `Commutative`.

use super::repr::Quaternion;
use crate::scalar::{Algebra, Conjugate, Field, Real, Ring};

impl<T: Real> Algebra for Quaternion<T> {
    const ZERO: Self = Self::new(T::ZERO, T::ZERO, T::ZERO, T::ZERO);
    const ONE: Self = Self::IDENTITY;
}

impl<T: Real> Ring for Quaternion<T> {}

impl<T: Real> Field for Quaternion<T> {
    #[inline]
    fn from_f64(value: f64) -> Self {
        Self::from_real(T::from_f64(value))
    }

    #[inline]
    fn recip(self) -> Self {
        self.inv()
    }
}

impl<T: Real> Conjugate for Quaternion<T> {
    type Real = T;

    #[inline]
    fn conj(self) -> Self {
        Self::new(self.w, -self.x, -self.y, -self.z)
    }

    #[inline]
    fn norm_sqr(self) -> T {
        self.dot(self)
    }

    #[inline]
    fn real_part(self) -> T {
        self.w
    }

    #[inline]
    fn max_abs(self) -> T {
        self.w.abs().max(self.x.abs()).max(self.y.abs().max(self.z.abs()))
    }

    #[inline]
    fn div_real(self, divisor: T) -> Self {
        self / divisor
    }

    #[inline]
    fn scale(self, factor: T) -> Self {
        self * factor
    }
}
