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

/// Division nach Smith ohne weitere Vorkehrung; läuft über, wenn die
/// Zwischenwerte es tun.
fn smith<T: Real>(lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
    if rhs.re.abs() >= rhs.im.abs() {
        let ratio = rhs.im / rhs.re;
        let denom = rhs.re + rhs.im * ratio;
        Complex::new((lhs.re + lhs.im * ratio) / denom, (lhs.im - lhs.re * ratio) / denom)
    } else {
        let ratio = rhs.re / rhs.im;
        let denom = rhs.re * ratio + rhs.im;
        Complex::new((lhs.re * ratio + lhs.im) / denom, (lhs.im * ratio - lhs.re) / denom)
    }
}

/// Division. Zähler und Nenner werden zuerst je auf die Größenordnung eins
/// gebracht, dann nach Smith geteilt; die Maßstäbe kommen am Ende dazu.
///
/// So laufen Zwischenwerte weder über noch unter: `z / z` ist auch für
/// Einträge nahe der größten oder kleinsten darstellbaren Zahl eins. Nur ein
/// Ergebnis, dessen Betrag der größten darstellbaren Zahl bis auf den Faktor
/// √2 nahe kommt, kann noch vorzeitig überlaufen. Division durch null
/// ergibt NaN-Anteile.
impl<T: Real> Div for Complex<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        let size = self.re.abs().max(self.im.abs());
        let scale = rhs.re.abs().max(rhs.im.abs());
        let usable = |value: T| value > T::ZERO && value.is_finite();
        if !(usable(size) && usable(scale)) {
            // Null, unendlich oder NaN: Hier gibt es nichts zu retten, die
            // einfache Formel liefert die üblichen Grenzwerte.
            return smith(self, rhs);
        }
        let quotient =
            smith(Self::new(self.re / size, self.im / size), Self::new(rhs.re / scale, rhs.im / scale));
        if scale >= T::ONE || size <= T::ONE {
            // Das Verhältnis der Maßstäbe läuft hier nicht vorzeitig über.
            let factor = size / scale;
            Self::new(quotient.re * factor, quotient.im * factor)
        } else {
            Self::new(quotient.re * size / scale, quotient.im * size / scale)
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
