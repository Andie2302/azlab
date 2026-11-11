//! Implementierungen der backend-freien Traits für die primitiven Typen.

use super::algebra::{Algebra, Commutative};
use super::conjugate::Conjugate;
use super::field::Field;
use super::real::Real;
use super::ring::Ring;

macro_rules! impl_ring {
    ($zero:literal, $one:literal => $($t:ty),*) => {
        $(
            impl Algebra for $t {
                const ZERO: Self = $zero;
                const ONE: Self = $one;
            }

            impl Ring for $t {}

            impl Commutative for $t {}
        )*
    };
}

impl_ring!(0, 1 => i8, i16, i32, i64, i128, isize);
impl_ring!(0.0, 1.0 => f32, f64);

macro_rules! impl_real {
    ($($t:ident),*) => {
        $(
            impl Field for $t {
                #[inline]
                fn from_f64(value: f64) -> Self {
                    value as $t
                }
            }

            impl Conjugate for $t {
                type Real = $t;

                #[inline]
                fn conj(self) -> Self {
                    self
                }
                #[inline]
                fn norm_sqr(self) -> Self {
                    self * self
                }
                #[inline]
                fn real_part(self) -> Self {
                    self
                }
                #[inline]
                fn scale(self, factor: Self) -> Self {
                    self * factor
                }
            }

            impl Real for $t {
                const PI: Self = core::$t::consts::PI;
                const E: Self = core::$t::consts::E;
                const EPSILON: Self = $t::EPSILON;
                const INFINITY: Self = $t::INFINITY;
                const NAN: Self = $t::NAN;

                #[inline]
                fn to_f64(self) -> f64 {
                    self as f64
                }
                #[inline]
                fn abs(self) -> Self {
                    $t::abs(self)
                }
                #[inline]
                fn min(self, other: Self) -> Self {
                    $t::min(self, other)
                }
                #[inline]
                fn max(self, other: Self) -> Self {
                    $t::max(self, other)
                }
                #[inline]
                fn copysign(self, sign: Self) -> Self {
                    $t::copysign(self, sign)
                }
                #[inline]
                fn is_finite(self) -> bool {
                    $t::is_finite(self)
                }
                #[inline]
                fn is_nan(self) -> bool {
                    $t::is_nan(self)
                }
            }
        )*
    };
}

impl_real!(f32, f64);
