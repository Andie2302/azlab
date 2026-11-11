//! Der Typ [`Complex`] und seine Grundfunktionen.

use crate::scalar::{Real, RealElementary};

/// Komplexe Zahl `re + im·i`.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Complex<T> {
    /// Realteil.
    pub re: T,
    /// Imaginärteil.
    pub im: T,
}

impl<T> Complex<T> {
    /// Erzeugt `re + im·i`.
    #[inline]
    pub const fn new(re: T, im: T) -> Self {
        Self { re, im }
    }
}

impl<T: Real> Complex<T> {
    /// Imaginäre Einheit `i`.
    pub const I: Self = Self::new(T::ZERO, T::ONE);

    /// Reelle Zahl als komplexe Zahl.
    #[inline]
    pub fn from_real(re: T) -> Self {
        Self::new(re, T::ZERO)
    }

    /// Multiplikation mit `i` (Drehung um 90°), ohne Rundungsfehler.
    #[inline]
    pub fn mul_i(self) -> Self {
        Self::new(-self.im, self.re)
    }

    /// Multiplikation mit `-i`, ohne Rundungsfehler.
    #[inline]
    pub fn mul_neg_i(self) -> Self {
        Self::new(self.im, -self.re)
    }
}

impl<T: RealElementary> Complex<T> {
    /// Betrag, überlaufsicher berechnet.
    #[inline]
    pub fn norm(self) -> T {
        self.re.hypot(self.im)
    }

    /// Argument (Winkel) im Bereich (−π, π].
    #[inline]
    pub fn arg(self) -> T {
        self.im.atan2(self.re)
    }

    /// Komplexe Zahl aus Betrag und Winkel.
    #[inline]
    pub fn from_polar(norm: T, arg: T) -> Self {
        let (sin, cos) = arg.sin_cos();
        Self::new(norm * cos, norm * sin)
    }

    /// Zerlegung in `(Betrag, Winkel)`.
    #[inline]
    pub fn to_polar(self) -> (T, T) {
        (self.norm(), self.arg())
    }
}
