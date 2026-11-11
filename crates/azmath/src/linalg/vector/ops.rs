//! Rechenoperatoren für [`Vector`].

use super::repr::Vector;
use crate::scalar::{Algebra, Field};
use core::ops::{Add, Div, DivAssign, Mul, MulAssign, Neg, Sub};

impl<T: Algebra, const N: usize> Add for Vector<T, N> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        self.zip_map(rhs, |a, b| a + b)
    }
}

impl<T: Algebra, const N: usize> Sub for Vector<T, N> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        self.zip_map(rhs, |a, b| a - b)
    }
}

impl<T: Algebra, const N: usize> Neg for Vector<T, N> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        self.map(|a| -a)
    }
}

/// Multiplikation jeder Komponente mit einem Skalar (von rechts).
impl<T: Algebra, const N: usize> Mul<T> for Vector<T, N> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        self.map(|a| a * rhs)
    }
}

/// Division jeder Komponente durch einen Skalar.
impl<T: Field, const N: usize> Div<T> for Vector<T, N> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self {
        self.map(|a| a / rhs)
    }
}

impl_op_assign!([T: Algebra, const N: usize] Vector<T, N>, AddAssign, add_assign, +);
impl_op_assign!([T: Algebra, const N: usize] Vector<T, N>, SubAssign, sub_assign, -);

impl<T: Algebra, const N: usize> MulAssign<T> for Vector<T, N> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

impl<T: Field, const N: usize> DivAssign<T> for Vector<T, N> {
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}

// `T * Vector<T, N>` lässt sich wegen der Kohärenzregeln nicht generisch
// schreiben, deshalb je primitivem Typ.
macro_rules! impl_left_scalar {
    ($($t:ty),*) => {
        $(
            impl<const N: usize> Mul<Vector<$t, N>> for $t {
                type Output = Vector<$t, N>;

                #[inline]
                fn mul(self, rhs: Vector<$t, N>) -> Vector<$t, N> {
                    rhs * self
                }
            }
        )*
    };
}

impl_left_scalar!(f32, f64);
