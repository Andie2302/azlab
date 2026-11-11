//! Rechenoperatoren für [`Interval`], alle mit Rundung nach außen.

use super::repr::Interval;
use crate::scalar::Real;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// Produkt zweier Grenzen mit der Intervall-Konvention `0 · ∞ = 0`.
#[inline]
fn bound_product<T: Real>(a: T, b: T) -> T {
    if a == T::ZERO || b == T::ZERO { T::ZERO } else { a * b }
}

impl<T: Real> Add for Interval<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::outward(self.lo + rhs.lo, self.hi + rhs.hi)
    }
}

impl<T: Real> Sub for Interval<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::outward(self.lo - rhs.hi, self.hi - rhs.lo)
    }
}

/// Exakt: Die Grenzen tauschen nur Platz und Vorzeichen.
impl<T: Real> Neg for Interval<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self { lo: -self.hi, hi: -self.lo }
    }
}

impl<T: Real> Mul for Interval<T> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let p1 = bound_product(self.lo, rhs.lo);
        let p2 = bound_product(self.lo, rhs.hi);
        let p3 = bound_product(self.hi, rhs.lo);
        let p4 = bound_product(self.hi, rhs.hi);
        Self::outward(p1.min(p2).min(p3.min(p4)), p1.max(p2).max(p3.max(p4)))
    }
}

/// Enthält der Divisor die Null, ist das Ergebnis die ganze Gerade.
impl<T: Real> Div for Interval<T> {
    type Output = Self;

    // Division als Multiplikation mit dem Kehrwert-Intervall.
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, rhs: Self) -> Self {
        if rhs.contains_zero() {
            return Self::ENTIRE;
        }
        let reciprocal = Self::outward(T::ONE / rhs.hi, T::ONE / rhs.lo);
        self * reciprocal
    }
}

impl_op_assign!([T: Real] Interval<T>, AddAssign, add_assign, +);
impl_op_assign!([T: Real] Interval<T>, SubAssign, sub_assign, -);
impl_op_assign!([T: Real] Interval<T>, MulAssign, mul_assign, *);
impl_op_assign!([T: Real] Interval<T>, DivAssign, div_assign, /);
