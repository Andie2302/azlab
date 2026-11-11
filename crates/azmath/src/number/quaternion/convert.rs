//! Umwandlungen von und nach [`Quaternion`].
//!
//! Die Zerlegung in zwei komplexe Zahlen ist `q = a + b·j` mit
//! `a = w + x·i` und `b = y + z·i`; sie passt zur Multiplikation der
//! Cayley-Dickson-Konstruktion.

use super::repr::Quaternion;
use crate::number::{CayleyDickson, Complex};
use crate::scalar::Real;

impl<T: Real> From<T> for Quaternion<T> {
    #[inline]
    fn from(w: T) -> Self {
        Self::from_real(w)
    }
}

impl<T: Real> From<Complex<T>> for Quaternion<T> {
    #[inline]
    fn from(z: Complex<T>) -> Self {
        Self::new(z.re, z.im, T::ZERO, T::ZERO)
    }
}

impl<T> From<CayleyDickson<CayleyDickson<T>>> for Quaternion<T> {
    #[inline]
    fn from(q: CayleyDickson<CayleyDickson<T>>) -> Self {
        Self::new(q.a.a, q.a.b, q.b.a, q.b.b)
    }
}

impl<T> From<Quaternion<T>> for CayleyDickson<CayleyDickson<T>> {
    #[inline]
    fn from(q: Quaternion<T>) -> Self {
        Self::new(CayleyDickson::new(q.w, q.x), CayleyDickson::new(q.y, q.z))
    }
}

impl<T> Quaternion<T> {
    /// Setzt `a + b·j` aus zwei komplexen Zahlen zusammen.
    #[inline]
    pub fn from_complex_pair(a: Complex<T>, b: Complex<T>) -> Self {
        Self::new(a.re, a.im, b.re, b.im)
    }

    /// Zerlegt in `(a, b)` mit `q = a + b·j`.
    #[inline]
    pub fn to_complex_pair(self) -> (Complex<T>, Complex<T>) {
        (Complex::new(self.w, self.x), Complex::new(self.y, self.z))
    }
}
