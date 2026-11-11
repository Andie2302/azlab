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
    /// Multiplikation jedes Anteils mit einer reellen Zahl.
    fn scale(self, factor: Self::Real) -> Self;
}
