//! Der Trait [`Conjugate`].

use super::algebra::Algebra;
use super::real::Real;

/// Konjugation und Norm.
///
/// Für reelle Typen ist die Konjugation die Identität; für komplexe Zahlen,
/// Quaternionen und allgemeine Cayley-Dickson-Algebren kehrt sie das
/// Vorzeichen der imaginären Anteile um.
pub trait Conjugate: Algebra {
    /// Der reelle Typ, in dem Normen gemessen werden.
    type Real: Real;

    /// Konjugierte Zahl.
    fn conj(self) -> Self;
    /// Quadrat der Norm, also `self * conj(self)` als reelle Zahl.
    fn norm_sqr(self) -> Self::Real;
    /// Realteil.
    fn real_part(self) -> Self::Real;
    /// Größter Betrag eines Anteils; für reelle Zahlen der Betrag.
    ///
    /// Ein Maß für die Größe wie die Norm, aber ohne Quadrieren und ohne
    /// Summieren: Es läuft weder über noch unter, solange die Zahl selbst
    /// darstellbar ist. Geeignet für Größenvergleiche (Pivotsuche) und zum
    /// Skalieren.
    fn max_abs(self) -> Self::Real;
    /// Division jedes Anteils durch eine reelle Zahl.
    ///
    /// Anders als `scale(1 / divisor)` bleibt das auch dann richtig, wenn der
    /// Kehrwert des Divisors nicht darstellbar ist.
    fn div_real(self, divisor: Self::Real) -> Self;
    /// Multiplikation jedes Anteils mit einer reellen Zahl.
    fn scale(self, factor: Self::Real) -> Self;
}
