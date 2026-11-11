//! Elementarfunktionen für [`Dual`] und [`DualN`] über die Kettenregel.
//!
//! An Stellen, an denen eine Funktion nicht differenzierbar ist (z. B.
//! `sqrt` bei 0), ergibt sich die Ableitung aus der Gleitkommarechnung und
//! kann unendlich oder NaN sein.

use super::repr::{Dual, DualN};
use super::tangent::Tangent;
use crate::scalar::{Elementary, RealElementary};

macro_rules! impl_dual_elementary {
    ([$($gen:tt)*] $ty:ty, $eps:ty) => {
        impl<$($gen)*> $ty
        where
            T: Elementary,
        {
            /// Kettenregel: neuer Wert, Ableitung mal `derivative`.
            #[inline]
            fn chain(self, value: T, derivative: T) -> Self {
                Self {
                    re: value,
                    eps: <$eps as Tangent<T>>::map(self.eps, |a| a * derivative),
                }
            }
        }

        impl<$($gen)*> Elementary for $ty
        where
            T: Elementary,
        {
            fn sqrt(self) -> Self {
                let value = self.re.sqrt();
                self.chain(value, (value + value).recip())
            }

            fn exp(self) -> Self {
                let value = self.re.exp();
                self.chain(value, value)
            }

            fn exp2(self) -> Self {
                let value = self.re.exp2();
                self.chain(value, value * T::from_f64(core::f64::consts::LN_2))
            }

            fn exp_m1(self) -> Self {
                self.chain(self.re.exp_m1(), self.re.exp())
            }

            fn ln(self) -> Self {
                self.chain(self.re.ln(), self.re.recip())
            }

            fn log2(self) -> Self {
                let scale = T::from_f64(core::f64::consts::LN_2);
                self.chain(self.re.log2(), (self.re * scale).recip())
            }

            fn log10(self) -> Self {
                let scale = T::from_f64(core::f64::consts::LN_10);
                self.chain(self.re.log10(), (self.re * scale).recip())
            }

            fn ln_1p(self) -> Self {
                self.chain(self.re.ln_1p(), (T::ONE + self.re).recip())
            }

            /// Ist der Exponent konstant, wird `ln` der Basis nicht gebraucht;
            /// so bleibt z. B. `x²` auch für negative `x` differenzierbar.
            fn powf(self, n: Self) -> Self {
                let value = self.re.powf(n.re);
                let by_base = n.re * self.re.powf(n.re - T::ONE);
                let by_exponent = value * self.re.ln();
                Self {
                    re: value,
                    eps: <$eps as Tangent<T>>::zip(self.eps, n.eps, |a, b| {
                        if b == T::ZERO {
                            a * by_base
                        } else {
                            a * by_base + b * by_exponent
                        }
                    }),
                }
            }

            fn sin(self) -> Self {
                self.chain(self.re.sin(), self.re.cos())
            }

            fn cos(self) -> Self {
                self.chain(self.re.cos(), -self.re.sin())
            }

            fn tan(self) -> Self {
                let value = self.re.tan();
                self.chain(value, T::ONE + value * value)
            }

            fn asin(self) -> Self {
                self.chain(self.re.asin(), (T::ONE - self.re * self.re).sqrt().recip())
            }

            fn acos(self) -> Self {
                self.chain(self.re.acos(), -(T::ONE - self.re * self.re).sqrt().recip())
            }

            fn atan(self) -> Self {
                self.chain(self.re.atan(), (T::ONE + self.re * self.re).recip())
            }

            fn sinh(self) -> Self {
                self.chain(self.re.sinh(), self.re.cosh())
            }

            fn cosh(self) -> Self {
                self.chain(self.re.cosh(), self.re.sinh())
            }

            fn tanh(self) -> Self {
                let value = self.re.tanh();
                self.chain(value, T::ONE - value * value)
            }

            fn asinh(self) -> Self {
                self.chain(self.re.asinh(), (self.re * self.re + T::ONE).sqrt().recip())
            }

            fn acosh(self) -> Self {
                self.chain(self.re.acosh(), (self.re * self.re - T::ONE).sqrt().recip())
            }

            fn atanh(self) -> Self {
                self.chain(self.re.atanh(), (T::ONE - self.re * self.re).recip())
            }
        }

        impl<$($gen)*> RealElementary for $ty
        where
            T: RealElementary,
        {
            fn cbrt(self) -> Self {
                let value = self.re.cbrt();
                self.chain(value, (T::from_f64(3.0) * value * value).recip())
            }

            fn atan2(self, x: Self) -> Self {
                let (y0, x0) = (self.re, x.re);
                let denom = x0 * x0 + y0 * y0;
                Self {
                    re: y0.atan2(x0),
                    eps: <$eps as Tangent<T>>::zip(self.eps, x.eps, |dy, dx| {
                        (x0 * dy - y0 * dx) / denom
                    }),
                }
            }

            fn hypot(self, other: Self) -> Self {
                let (a0, b0) = (self.re, other.re);
                let value = a0.hypot(b0);
                Self {
                    re: value,
                    eps: <$eps as Tangent<T>>::zip(self.eps, other.eps, |da, db| {
                        (a0 * da + b0 * db) / value
                    }),
                }
            }

            /// Stückweise konstant: Ableitung null.
            fn floor(self) -> Self {
                self.chain(self.re.floor(), T::ZERO)
            }

            /// Stückweise konstant: Ableitung null.
            fn ceil(self) -> Self {
                self.chain(self.re.ceil(), T::ZERO)
            }

            /// Stückweise konstant: Ableitung null.
            fn round(self) -> Self {
                self.chain(self.re.round(), T::ZERO)
            }

            /// Stückweise konstant: Ableitung null.
            fn trunc(self) -> Self {
                self.chain(self.re.trunc(), T::ZERO)
            }

            fn mul_add(self, a: Self, b: Self) -> Self {
                let (x, y) = (self.re, a.re);
                let product = <$eps as Tangent<T>>::zip(self.eps, a.eps, |dx, dy| x * dy + dx * y);
                Self {
                    re: x.mul_add(y, b.re),
                    eps: <$eps as Tangent<T>>::zip(product, b.eps, |p, db| p + db),
                }
            }
        }
    };
}

impl_dual_elementary!([T] Dual<T>, T);
impl_dual_elementary!([T, const N: usize] DualN<T, N>, [T; N]);
