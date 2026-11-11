//! Elementarfunktionen für [`Complex`].
//!
//! Es gelten die Hauptzweige wie in C99: `ln` und `sqrt` haben ihren Schnitt
//! entlang der negativen reellen Achse, das Argument liegt in (−π, π].
//!
//! Die Umkehrfunktionen (`asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`)
//! benutzen die Lehrbuchformeln über `ln` und `sqrt`. Sie sind korrekt, aber
//! nahe den Verzweigungsschnitten und für sehr große Argumente nicht auf
//! letzte Stelle genau.

use super::repr::Complex;
use crate::scalar::{Algebra, Elementary, RealElementary};

impl<T: RealElementary> Complex<T> {
    fn half() -> T {
        T::from_f64(0.5)
    }

    /// Stabile Form von `tanh` für die Zerlegung `x + y·i`.
    fn tanh_parts(x: T, y: T) -> Self {
        let two_x = x + x;
        let two_y = y + y;
        // Ab hier ist tanh(x) im Rahmen der Genauigkeit ±1 und cosh(2x)
        // würde unnötig groß.
        if two_x.abs() > T::from_f64(40.0) {
            Self::new(T::ONE.copysign(x), T::ZERO)
        } else {
            let denom = two_x.cosh() + two_y.cos();
            Self::new(two_x.sinh() / denom, two_y.sin() / denom)
        }
    }
}

impl<T: RealElementary> Elementary for Complex<T> {
    fn sqrt(self) -> Self {
        if self.re == T::ZERO && self.im == T::ZERO {
            return Self::new(T::ZERO, self.im);
        }
        let norm = self.norm();
        if self.re >= T::ZERO {
            let t = ((norm + self.re) * Self::half()).sqrt();
            Self::new(t, self.im / (t + t))
        } else {
            let t = ((norm - self.re) * Self::half()).sqrt();
            Self::new(self.im.abs() / (t + t), t.copysign(self.im))
        }
    }

    fn exp(self) -> Self {
        let (sin, cos) = self.im.sin_cos();
        let radius = self.re.exp();
        Self::new(radius * cos, radius * sin)
    }

    fn exp2(self) -> Self {
        (self * T::from_f64(core::f64::consts::LN_2)).exp()
    }

    fn exp_m1(self) -> Self {
        // e^a·cos b − 1 = expm1(a)·cos b + (cos b − 1), und
        // cos b − 1 = −2·sin²(b/2) vermeidet die Auslöschung.
        let (sin, cos) = self.im.sin_cos();
        let sin_half = (self.im * Self::half()).sin();
        let cos_m1 = -(sin_half * sin_half) - sin_half * sin_half;
        Self::new(self.re.exp_m1() * cos + cos_m1, self.re.exp() * sin)
    }

    fn ln(self) -> Self {
        Self::new(self.norm().ln(), self.arg())
    }

    fn log2(self) -> Self {
        self.ln() / T::from_f64(core::f64::consts::LN_2)
    }

    fn log10(self) -> Self {
        self.ln() / T::from_f64(core::f64::consts::LN_10)
    }

    fn ln_1p(self) -> Self {
        // |1 + z|² = 1 + a·(2 + a) + b²
        let two = T::ONE + T::ONE;
        let re = (self.re * (two + self.re) + self.im * self.im).ln_1p() * Self::half();
        Self::new(re, self.im.atan2(T::ONE + self.re))
    }

    fn powf(self, n: Self) -> Self {
        if self == Self::ZERO {
            return if n == Self::ZERO { Self::ONE } else { Self::ZERO };
        }
        (n * self.ln()).exp()
    }

    fn sin(self) -> Self {
        let (sin, cos) = self.re.sin_cos();
        Self::new(sin * self.im.cosh(), cos * self.im.sinh())
    }

    fn cos(self) -> Self {
        let (sin, cos) = self.re.sin_cos();
        Self::new(cos * self.im.cosh(), -(sin * self.im.sinh()))
    }

    fn tan(self) -> Self {
        // tan z = −i·tanh(i·z)
        Self::tanh_parts(-self.im, self.re).mul_neg_i()
    }

    fn asin(self) -> Self {
        // asin z = −i·asinh(i·z)
        self.mul_i().asinh().mul_neg_i()
    }

    fn acos(self) -> Self {
        Self::from_real(T::PI * Self::half()) - self.asin()
    }

    fn atan(self) -> Self {
        // atan z = −i·atanh(i·z)
        self.mul_i().atanh().mul_neg_i()
    }

    fn sinh(self) -> Self {
        let (sin, cos) = self.im.sin_cos();
        Self::new(self.re.sinh() * cos, self.re.cosh() * sin)
    }

    fn cosh(self) -> Self {
        let (sin, cos) = self.im.sin_cos();
        Self::new(self.re.cosh() * cos, self.re.sinh() * sin)
    }

    fn tanh(self) -> Self {
        Self::tanh_parts(self.re, self.im)
    }

    fn asinh(self) -> Self {
        (self + (self * self + Self::ONE).sqrt()).ln()
    }

    fn acosh(self) -> Self {
        (self + (self + Self::ONE).sqrt() * (self - Self::ONE).sqrt()).ln()
    }

    fn atanh(self) -> Self {
        ((Self::ONE + self).ln() - (Self::ONE - self).ln()) * Self::half()
    }
}
