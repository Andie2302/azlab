//! Rechenoperatoren für [`Quaternion`].

use super::repr::Quaternion;
use crate::scalar::Real;
use core::ops::{Add, Div, Mul, Neg, Sub};

impl<T: Real> Add for Quaternion<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.w + rhs.w, self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: Real> Sub for Quaternion<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.w - rhs.w, self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: Real> Neg for Quaternion<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.w, -self.x, -self.y, -self.z)
    }
}

/// Hamilton-Produkt; nicht kommutativ.
impl<T: Real> Mul for Quaternion<T> {
    type Output = Self;

    // Das Hamilton-Produkt besteht aus Summen und Differenzen von Produkten.
    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let (a, b) = (self, rhs);
        Self::new(
            a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
        )
    }
}

/// Rechtsdivision `self · rhs⁻¹`.
impl<T: Real> Div for Quaternion<T> {
    type Output = Self;

    // Division ist als Multiplikation mit dem Inversen definiert.
    #[allow(clippy::suspicious_arithmetic_impl)]
    #[inline]
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}

impl_op_assign!([T: Real] Quaternion<T>, AddAssign, add_assign, +);
impl_op_assign!([T: Real] Quaternion<T>, SubAssign, sub_assign, -);
impl_op_assign!([T: Real] Quaternion<T>, MulAssign, mul_assign, *);
impl_op_assign!([T: Real] Quaternion<T>, DivAssign, div_assign, /);

/// Skalierung mit einer reellen Zahl.
impl<T: Real> Mul<T> for Quaternion<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        Self::new(self.w * rhs, self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

/// Division durch eine reelle Zahl.
impl<T: Real> Div<T> for Quaternion<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self {
        Self::new(self.w / rhs, self.x / rhs, self.y / rhs, self.z / rhs)
    }
}

// `T * Quaternion<T>` lässt sich wegen der Kohärenzregeln nicht generisch
// schreiben, deshalb je primitivem Typ.
macro_rules! impl_left_scalar {
    ($($t:ty),*) => {
        $(
            impl Mul<Quaternion<$t>> for $t {
                type Output = Quaternion<$t>;

                #[inline]
                fn mul(self, rhs: Quaternion<$t>) -> Quaternion<$t> {
                    rhs * self
                }
            }
        )*
    };
}

impl_left_scalar!(f32, f64);
