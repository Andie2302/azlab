//! Umwandlungen zwischen [`CayleyDickson`] und [`Complex`].

use super::repr::CayleyDickson;
use crate::number::Complex;

impl<T> From<Complex<T>> for CayleyDickson<T> {
    #[inline]
    fn from(z: Complex<T>) -> Self {
        Self::new(z.re, z.im)
    }
}

impl<T> From<CayleyDickson<T>> for Complex<T> {
    #[inline]
    fn from(z: CayleyDickson<T>) -> Self {
        Self::new(z.a, z.b)
    }
}
