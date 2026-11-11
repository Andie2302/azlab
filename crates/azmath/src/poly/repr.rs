//! Der Typ [`Polynomial`].

use super::eval::{horner, horner_with_derivative};
use crate::scalar::Algebra;

/// Polynom mit `N` Koeffizienten (Grad höchstens `N − 1`), aufsteigend:
/// `c[0] + c[1]·x + … + c[N−1]·x^(N−1)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Polynomial<T, const N: usize>(pub [T; N]);

impl<T, const N: usize> Polynomial<T, N> {
    /// Erzeugt ein Polynom aus seinen Koeffizienten (aufsteigend).
    #[inline]
    pub const fn new(coefficients: [T; N]) -> Self {
        Self(coefficients)
    }

    /// Die Koeffizienten, aufsteigend.
    #[inline]
    pub fn coefficients(&self) -> &[T] {
        &self.0
    }

    /// Wandelt jeden Koeffizienten um, z. B. in Dualzahlen, komplexe Zahlen
    /// oder Intervalle, um das Polynom dort auszuwerten.
    #[inline]
    pub fn map<U>(self, f: impl FnMut(T) -> U) -> Polynomial<U, N> {
        Polynomial(self.0.map(f))
    }
}

impl<T: Algebra, const N: usize> Polynomial<T, N> {
    /// Das Nullpolynom.
    pub const ZERO: Self = Self([T::ZERO; N]);

    /// Wert an der Stelle `x` (Horner-Schema).
    #[inline]
    pub fn eval(&self, x: T) -> T {
        horner(&self.0, x)
    }

    /// Wert und Ableitung an der Stelle `x`: `(Wert, Ableitung)`.
    #[inline]
    pub fn eval_with_derivative(&self, x: T) -> (T, T) {
        horner_with_derivative(&self.0, x)
    }

    /// Ableitung als Polynom gleicher Größe; der höchste Koeffizient ist
    /// dann null.
    pub fn derivative(&self) -> Self {
        let mut out = [T::ZERO; N];
        let mut exponent = T::ZERO;
        for (slot, &c) in out.iter_mut().zip(self.0.iter().skip(1)) {
            exponent += T::ONE;
            *slot = c * exponent;
        }
        Self(out)
    }

    /// Grad: Index des höchsten Koeffizienten ungleich null. `None` für das
    /// Nullpolynom.
    #[inline]
    pub fn degree(&self) -> Option<usize> {
        self.0.iter().rposition(|c| *c != T::ZERO)
    }
}
