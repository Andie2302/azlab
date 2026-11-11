//! Kommutator und Exponentialfunktion quadratischer Matrizen.

use super::repr::Matrix;
use crate::scalar::{Algebra, Conjugate, Field, RealElementary};

/// Anzahl der Glieder der Taylor-Reihe nach dem Skalieren.
const TERMS: usize = 20;

/// Höchstzahl der Halbierungen beim Skalieren.
const MAX_SQUARINGS: u32 = 64;

impl<T: Algebra, const N: usize> Matrix<T, N, N> {
    /// Kommutator `self·other − other·self`: die Lie-Klammer von Matrizen.
    #[inline]
    pub fn commutator(self, other: Self) -> Self {
        self * other - other * self
    }
}

impl<T: Field + Conjugate, const N: usize> Matrix<T, N, N>
where
    T::Real: RealElementary,
{
    /// Matrix-Exponential `e^A`.
    ///
    /// Skalieren und Quadrieren: Die Matrix wird halbiert, bis ihre Norm
    /// klein ist, dann die Taylor-Reihe ausgewertet und das Ergebnis
    /// entsprechend oft quadriert. Gedacht für kleine Matrizen.
    pub fn exp(self) -> Self {
        let half = <T::Real as Field>::from_f64(0.5);
        let mut scaled = self;
        let mut squarings = 0;
        while scaled.norm() > half && squarings < MAX_SQUARINGS {
            scaled = scaled.map(|x| x.scale(half));
            squarings += 1;
        }

        let mut sum = Self::identity();
        let mut term = Self::identity();
        let mut index = T::ZERO;
        for _ in 0..TERMS {
            index += T::ONE;
            term = (term * scaled) / index;
            sum += term;
        }
        for _ in 0..squarings {
            sum = sum * sum;
        }
        sum
    }
}
