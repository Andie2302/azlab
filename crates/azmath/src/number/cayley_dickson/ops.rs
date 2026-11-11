//! Rechenoperatoren für [`CayleyDickson`].

use super::repr::CayleyDickson;
use crate::scalar::Conjugate;
use core::ops::{Add, Div, Mul, Neg, Sub};

impl<A: Conjugate> Add for CayleyDickson<A> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.a + rhs.a, self.b + rhs.b)
    }
}

impl<A: Conjugate> Sub for CayleyDickson<A> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.a - rhs.a, self.b - rhs.b)
    }
}

impl<A: Conjugate> Neg for CayleyDickson<A> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.a, -self.b)
    }
}

/// `(a, b)·(c, d) = (a·c − d*·b, d·a + b·c*)`
impl<A: Conjugate> Mul for CayleyDickson<A> {
    type Output = Self;

    // Die Formel besteht aus Summe und Differenz von Produkten.
    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.a * rhs.a - rhs.b.conj() * self.b,
            rhs.b * self.a + self.b * rhs.a.conj(),
        )
    }
}

/// Rechtsdivision `self · rhs⁻¹`.
impl<A: Conjugate> Div for CayleyDickson<A> {
    type Output = Self;

    // Division ist als Multiplikation mit dem Inversen definiert.
    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}

impl_op_assign!([A: Conjugate] CayleyDickson<A>, AddAssign, add_assign, +);
impl_op_assign!([A: Conjugate] CayleyDickson<A>, SubAssign, sub_assign, -);
impl_op_assign!([A: Conjugate] CayleyDickson<A>, MulAssign, mul_assign, *);
impl_op_assign!([A: Conjugate] CayleyDickson<A>, DivAssign, div_assign, /);
