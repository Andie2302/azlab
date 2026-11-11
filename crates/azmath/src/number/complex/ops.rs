//! Rechenoperatoren für [`Complex`].

use super::repr::Complex;
use crate::scalar::Real;
use core::ops::{Add, Div, Mul, Neg, Sub};

impl<T: Real> Add for Complex<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl<T: Real> Sub for Complex<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl<T: Real> Neg for Complex<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

impl<T: Real> Mul for Complex<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

/// Division nach Smith: vermeidet Überlauf, solange das Ergebnis darstellbar
/// ist. Division durch null ergibt NaN-Anteile.
impl<T: Real> Div for Complex<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        if rhs.re.abs() >= rhs.im.abs() {
            let ratio = rhs.im / rhs.re;
            let denom = rhs.re + rhs.im * ratio;
            Self::new(
                (self.re + self.im * ratio) / denom,
                (self.im - self.re * ratio) / denom,
            )
        } else {
            let ratio = rhs.re / rhs.im;
            let denom = rhs.re * ratio + rhs.im;
            Self::new(
                (self.re * ratio + self.im) / denom,
                (self.im * ratio - self.re) / denom,
            )
        }
    }
}

impl_op_assign!([T: Real] Complex<T>, AddAssign, add_assign, +);
impl_op_assign!([T: Real] Complex<T>, SubAssign, sub_assign, -);
impl_op_assign!([T: Real] Complex<T>, MulAssign, mul_assign, *);
impl_op_assign!([T: Real] Complex<T>, DivAssign, div_assign, /);

/// Skalierung mit einer reellen Zahl.
impl<T: Real> Mul<T> for Complex<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

/// Division durch eine reelle Zahl.
impl<T: Real> Div<T> for Complex<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self {
        Self::new(self.re / rhs, self.im / rhs)
    }
}

// `T * Complex<T>` lässt sich wegen der Kohärenzregeln nicht generisch
// schreiben, deshalb je primitivem Typ.
macro_rules! impl_left_scalar {
    ($($t:ty),*) => {
        $(
            impl Mul<Complex<$t>> for $t {
                type Output = Complex<$t>;

                #[inline]
                fn mul(self, rhs: Complex<$t>) -> Complex<$t> {
                    rhs * self
                }
            }
        )*
    };
}

impl_left_scalar!(f32, f64);
