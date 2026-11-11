//! Skalar-Traits für [`Dual`] und [`DualN`].

use super::repr::{Dual, DualN};
use super::tangent::Tangent;
use crate::scalar::{Algebra, Commutative, Conjugate, Field, Real, Ring};

macro_rules! impl_dual_traits {
    ([$($gen:tt)*] $ty:ty, $eps:ty, $zero:expr) => {
        impl<$($gen)*> Algebra for $ty
        where
            T: Algebra,
        {
            const ZERO: Self = Self { re: T::ZERO, eps: $zero };
            const ONE: Self = Self { re: T::ONE, eps: $zero };
        }

        impl<$($gen)*> Ring for $ty where T: Ring {}

        impl<$($gen)*> Commutative for $ty where T: Commutative {}

        impl<$($gen)*> Field for $ty
        where
            T: Field,
        {
            #[inline]
            fn from_f64(value: f64) -> Self {
                Self { re: T::from_f64(value), eps: $zero }
            }
        }

        /// Dualzahlen über reellen Zahlen sind selbst reell: Die Konjugation
        /// hat keine Wirkung, Normen bleiben Dualzahlen und tragen die
        /// Ableitung weiter.
        impl<$($gen)*> Conjugate for $ty
        where
            T: Real,
        {
            type Real = Self;

            #[inline]
            fn conj(self) -> Self {
                self
            }

            #[inline]
            fn norm_sqr(self) -> Self {
                self * self
            }

            #[inline]
            fn real_part(self) -> Self {
                self
            }

            #[inline]
            fn max_abs(self) -> Self {
                <Self as Real>::abs(self)
            }

            #[inline]
            fn div_real(self, divisor: Self) -> Self {
                self / divisor
            }

            #[inline]
            fn scale(self, factor: Self) -> Self {
                self * factor
            }
        }

        /// Alle Abfragen und Vergleiche beziehen sich nur auf den Wert.
        impl<$($gen)*> Real for $ty
        where
            T: Real,
        {
            const PI: Self = Self { re: T::PI, eps: $zero };
            const E: Self = Self { re: T::E, eps: $zero };
            const EPSILON: Self = Self { re: T::EPSILON, eps: $zero };
            const INFINITY: Self = Self { re: T::INFINITY, eps: $zero };
            const NEG_INFINITY: Self = Self { re: T::NEG_INFINITY, eps: $zero };
            const NAN: Self = Self { re: T::NAN, eps: $zero };

            #[inline]
            fn to_f64(self) -> f64 {
                self.re.to_f64()
            }

            #[inline]
            fn abs(self) -> Self {
                if self.re < T::ZERO { -self } else { self }
            }

            #[inline]
            fn min(self, other: Self) -> Self {
                if other.re < self.re || self.re.is_nan() { other } else { self }
            }

            #[inline]
            fn max(self, other: Self) -> Self {
                if other.re > self.re || self.re.is_nan() { other } else { self }
            }

            #[inline]
            fn copysign(self, sign: Self) -> Self {
                let re = self.re.copysign(sign.re);
                let eps = if re == self.re {
                    self.eps
                } else {
                    <$eps as Tangent<T>>::map(self.eps, |a| -a)
                };
                Self { re, eps }
            }

            #[inline]
            fn next_up(self) -> Self {
                Self { re: self.re.next_up(), eps: self.eps }
            }

            #[inline]
            fn next_down(self) -> Self {
                Self { re: self.re.next_down(), eps: self.eps }
            }

            #[inline]
            fn is_finite(self) -> bool {
                self.re.is_finite()
            }

            #[inline]
            fn is_nan(self) -> bool {
                self.re.is_nan()
            }
        }
    };
}

impl_dual_traits!([T] Dual<T>, T, T::ZERO);
impl_dual_traits!([T, const N: usize] DualN<T, N>, [T; N], [T::ZERO; N]);
