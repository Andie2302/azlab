//! Vergleiche und Rechenoperatoren für [`Dual`] und [`DualN`].

use super::repr::{Dual, DualN};
use super::tangent::Tangent;
use crate::scalar::{Algebra, Field};
use core::cmp::Ordering;
use core::ops::{Add, Div, Mul, Neg, Sub};

macro_rules! impl_dual_ops {
    ([$($gen:tt)*] $ty:ty, $eps:ty) => {
        /// Vergleicht nur den Wert, nicht die Ableitung.
        impl<$($gen)*> PartialEq for $ty
        where
            T: PartialEq,
        {
            #[inline]
            fn eq(&self, other: &Self) -> bool {
                self.re == other.re
            }
        }

        /// Ordnet nur nach dem Wert, nicht nach der Ableitung.
        impl<$($gen)*> PartialOrd for $ty
        where
            T: PartialOrd,
        {
            #[inline]
            fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
                self.re.partial_cmp(&other.re)
            }
        }

        impl<$($gen)*> Add for $ty
        where
            T: Algebra,
        {
            type Output = Self;

            #[inline]
            fn add(self, rhs: Self) -> Self {
                Self {
                    re: self.re + rhs.re,
                    eps: <$eps as Tangent<T>>::zip(self.eps, rhs.eps, |a, b| a + b),
                }
            }
        }

        impl<$($gen)*> Sub for $ty
        where
            T: Algebra,
        {
            type Output = Self;

            #[inline]
            fn sub(self, rhs: Self) -> Self {
                Self {
                    re: self.re - rhs.re,
                    eps: <$eps as Tangent<T>>::zip(self.eps, rhs.eps, |a, b| a - b),
                }
            }
        }

        impl<$($gen)*> Neg for $ty
        where
            T: Algebra,
        {
            type Output = Self;

            #[inline]
            fn neg(self) -> Self {
                Self {
                    re: -self.re,
                    eps: <$eps as Tangent<T>>::map(self.eps, |a| -a),
                }
            }
        }

        /// Produktregel; die Reihenfolge der Faktoren bleibt erhalten.
        impl<$($gen)*> Mul for $ty
        where
            T: Algebra,
        {
            type Output = Self;

            // Produktregel: die Summe ist hier Absicht.
            #[allow(clippy::suspicious_arithmetic_impl)]
            #[inline]
            fn mul(self, rhs: Self) -> Self {
                let (x, y) = (self.re, rhs.re);
                Self {
                    re: x * y,
                    eps: <$eps as Tangent<T>>::zip(self.eps, rhs.eps, |a, b| x * b + a * y),
                }
            }
        }

        /// Quotientenregel in der Form `(a' − q·b') / b` mit `q = a / b`.
        impl<$($gen)*> Div for $ty
        where
            T: Field,
        {
            type Output = Self;

            // Quotientenregel: Differenz und Produkt sind hier Absicht.
            #[allow(clippy::suspicious_arithmetic_impl)]
            #[inline]
            fn div(self, rhs: Self) -> Self {
                let y = rhs.re;
                let q = self.re / y;
                Self {
                    re: q,
                    eps: <$eps as Tangent<T>>::zip(self.eps, rhs.eps, |a, b| (a - q * b) / y),
                }
            }
        }

        impl_op_assign!([$($gen)*] $ty, AddAssign, add_assign, +, where T: Algebra);
        impl_op_assign!([$($gen)*] $ty, SubAssign, sub_assign, -, where T: Algebra);
        impl_op_assign!([$($gen)*] $ty, MulAssign, mul_assign, *, where T: Algebra);
        impl_op_assign!([$($gen)*] $ty, DivAssign, div_assign, /, where T: Field);
    };
}

impl_dual_ops!([T] Dual<T>, T);
impl_dual_ops!([T, const N: usize] DualN<T, N>, [T; N]);
