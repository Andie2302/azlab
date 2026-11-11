//! Umwandlungen von und nach [`Complex`].

use super::repr::Complex;
use crate::scalar::Real;

impl<T: Real> From<T> for Complex<T> {
    #[inline]
    fn from(re: T) -> Self {
        Self::from_real(re)
    }
}

impl<T> From<(T, T)> for Complex<T> {
    #[inline]
    fn from((re, im): (T, T)) -> Self {
        Self::new(re, im)
    }
}

impl<T> From<Complex<T>> for (T, T) {
    #[inline]
    fn from(z: Complex<T>) -> Self {
        (z.re, z.im)
    }
}
