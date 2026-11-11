//! Die Traits [`Elementary`] und [`RealElementary`].
//!
//! Für `f32`/`f64` stammen die Implementierungen aus dem per Feature
//! gewählten Backend (`libm` oder `std`). Eigene Zahlentypen (z. B.
//! Dualzahlen) implementieren die Traits selbst und sind damit in allen
//! Algorithmen einsetzbar, die sie fordern.

use super::field::Field;
use super::real::Real;

/// Elementarfunktionen, die keine Ordnung voraussetzen.
///
/// Dadurch auch für komplexe Zahlen implementierbar.
pub trait Elementary: Field {
    /// Quadratwurzel.
    fn sqrt(self) -> Self;
    /// Exponentialfunktion `e^self`.
    fn exp(self) -> Self;
    /// `2^self`.
    fn exp2(self) -> Self;
    /// `e^self - 1`, genau auch für kleine Argumente.
    fn exp_m1(self) -> Self;
    /// Natürlicher Logarithmus.
    fn ln(self) -> Self;
    /// Logarithmus zur Basis 2.
    fn log2(self) -> Self;
    /// Logarithmus zur Basis 10.
    fn log10(self) -> Self;
    /// `ln(1 + self)`, genau auch für kleine Argumente.
    fn ln_1p(self) -> Self;
    /// Potenz mit Exponenten desselben Typs.
    fn powf(self, n: Self) -> Self;
    /// Sinus (Bogenmaß).
    fn sin(self) -> Self;
    /// Kosinus (Bogenmaß).
    fn cos(self) -> Self;
    /// Tangens (Bogenmaß).
    fn tan(self) -> Self;
    /// Arkussinus.
    fn asin(self) -> Self;
    /// Arkuskosinus.
    fn acos(self) -> Self;
    /// Arkustangens.
    fn atan(self) -> Self;
    /// Sinus hyperbolicus.
    fn sinh(self) -> Self;
    /// Kosinus hyperbolicus.
    fn cosh(self) -> Self;
    /// Tangens hyperbolicus.
    fn tanh(self) -> Self;
    /// Areasinus hyperbolicus.
    fn asinh(self) -> Self;
    /// Areakosinus hyperbolicus.
    fn acosh(self) -> Self;
    /// Areatangens hyperbolicus.
    fn atanh(self) -> Self;

    /// Sinus und Kosinus zugleich als `(sin, cos)`.
    #[inline]
    fn sin_cos(self) -> (Self, Self) {
        (self.sin(), self.cos())
    }
}

/// Elementarfunktionen, die nur für reelle Zahlen sinnvoll sind.
pub trait RealElementary: Real + Elementary {
    /// Kubikwurzel.
    fn cbrt(self) -> Self;
    /// Arkustangens von `self / x` mit korrektem Quadranten (`self` ist y).
    fn atan2(self, x: Self) -> Self;
    /// `sqrt(self² + other²)` ohne vermeidbaren Überlauf.
    fn hypot(self, other: Self) -> Self;
    /// Größte ganze Zahl kleiner oder gleich `self`.
    fn floor(self) -> Self;
    /// Kleinste ganze Zahl größer oder gleich `self`.
    fn ceil(self) -> Self;
    /// Rundet zur nächsten ganzen Zahl, bei ,5 von der Null weg.
    fn round(self) -> Self;
    /// Schneidet den Nachkommaanteil ab.
    fn trunc(self) -> Self;
    /// `self * a + b` mit nur einer Rundung.
    fn mul_add(self, a: Self, b: Self) -> Self;
}
