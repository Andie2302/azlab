//! Skalar-Traits für [`CayleyDickson`].
//!
//! Die Schranken sorgen dafür, dass jede Stufe genau die Eigenschaften
//! zusagt, die sie mathematisch hat.

use super::repr::CayleyDickson;
use crate::scalar::{Algebra, Commutative, Conjugate, Field, Real, Ring};

impl<A: Conjugate> Algebra for CayleyDickson<A> {
    const ZERO: Self = Self::new(A::ZERO, A::ZERO);
    const ONE: Self = Self::new(A::ONE, A::ZERO);
}

impl<A: Conjugate> Conjugate for CayleyDickson<A> {
    type Real = A::Real;

    #[inline]
    fn conj(self) -> Self {
        Self::new(self.a.conj(), -self.b)
    }

    #[inline]
    fn norm_sqr(self) -> A::Real {
        self.a.norm_sqr() + self.b.norm_sqr()
    }

    #[inline]
    fn real_part(self) -> A::Real {
        self.a.real_part()
    }

    #[inline]
    fn scale(self, factor: A::Real) -> Self {
        Self::new(self.a.scale(factor), self.b.scale(factor))
    }
}

/// Die Verdopplung ist genau dann assoziativ, wenn die Ausgangsalgebra
/// assoziativ und kommutativ ist.
impl<A: Ring + Commutative + Conjugate> Ring for CayleyDickson<A> {}

impl<A: Ring + Commutative + Conjugate> Field for CayleyDickson<A> {
    #[inline]
    fn from_f64(value: f64) -> Self {
        Self::ONE.scale(<A::Real as Field>::from_f64(value))
    }
}

/// Die Verdopplung ist genau dann kommutativ, wenn die Ausgangsalgebra
/// reell ist (Konjugation ohne Wirkung).
impl<A: Real + Conjugate> Commutative for CayleyDickson<A> {}
