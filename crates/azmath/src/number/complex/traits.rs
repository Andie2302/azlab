//! Skalar-Traits für [`Complex`].

use super::repr::Complex;
use crate::scalar::{Algebra, Commutative, Conjugate, Field, Real, Ring};

impl<T: Real> Algebra for Complex<T> {
    const ZERO: Self = Self::new(T::ZERO, T::ZERO);
    const ONE: Self = Self::new(T::ONE, T::ZERO);
}

impl<T: Real> Ring for Complex<T> {}

impl<T: Real> Commutative for Complex<T> {}

impl<T: Real> Field for Complex<T> {
    #[inline]
    fn from_f64(value: f64) -> Self {
        Self::from_real(T::from_f64(value))
    }
}

impl<T: Real> Conjugate for Complex<T> {
    type Real = T;

    #[inline]
    fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    #[inline]
    fn norm_sqr(self) -> T {
        self.re * self.re + self.im * self.im
    }

    #[inline]
    fn real_part(self) -> T {
        self.re
    }

    #[inline]
    fn scale(self, factor: T) -> Self {
        self * factor
    }
}
