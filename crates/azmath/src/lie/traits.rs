//! Die Traits [`LieAlgebra`] und [`LieGroup`].

use crate::scalar::Real;
use core::ops::{Add, Neg, Sub};

/// Lie-Algebra: Vektorraum mit Lie-Klammer.
pub trait LieAlgebra: Copy + Add<Output = Self> + Sub<Output = Self> + Neg<Output = Self> {
    /// Der reelle Skalartyp.
    type Scalar: Real;

    /// Nullvektor.
    const ZERO: Self;

    /// Multiplikation mit einem Skalar.
    fn scale(self, factor: Self::Scalar) -> Self;

    /// Lie-Klammer `[self, other]`.
    fn bracket(self, other: Self) -> Self;
}

/// Lie-Gruppe mit Exponential- und Logarithmusabbildung.
pub trait LieGroup: Copy {
    /// Die zugehörige Lie-Algebra (Tangentialraum am neutralen Element).
    type Algebra: LieAlgebra;

    /// Neutrales Element.
    fn identity() -> Self;

    /// Verknüpfung: erst `other`, dann `self`.
    fn compose(self, other: Self) -> Self;

    /// Inverses Element.
    fn inverse(self) -> Self;

    /// Exponentialabbildung: Algebra → Gruppe.
    fn exp(tangent: Self::Algebra) -> Self;

    /// Logarithmus: Gruppe → Algebra (Hauptwert).
    fn log(self) -> Self::Algebra;

    /// Adjungierte Darstellung `Ad_g(ξ)`: Es gilt
    /// `g · exp(ξ) · g⁻¹ = exp(Ad_g(ξ))`.
    fn adjoint(self, tangent: Self::Algebra) -> Self::Algebra;
}
