//! Rechenoperatoren für [`Polynomial`].

use super::repr::Polynomial;
use crate::scalar::Algebra;
use core::ops::{Add, Mul, Neg, Sub};

impl<T: Algebra, const N: usize> Add for Polynomial<T, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut out = self.0;
        for (slot, b) in out.iter_mut().zip(rhs.0) {
            *slot += b;
        }
        Self(out)
    }
}

impl<T: Algebra, const N: usize> Sub for Polynomial<T, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut out = self.0;
        for (slot, b) in out.iter_mut().zip(rhs.0) {
            *slot -= b;
        }
        Self(out)
    }
}

impl<T: Algebra, const N: usize> Neg for Polynomial<T, N> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        self.map(|c| -c)
    }
}

/// Multiplikation jedes Koeffizienten mit einem Skalar.
impl<T: Algebra, const N: usize> Mul<T> for Polynomial<T, N> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        self.map(|c| c * rhs)
    }
}

impl_op_assign!([T: Algebra, const N: usize] Polynomial<T, N>, AddAssign, add_assign, +);
impl_op_assign!([T: Algebra, const N: usize] Polynomial<T, N>, SubAssign, sub_assign, -);
