//! Produkte und Normen für [`Vector`].

use super::repr::Vector;
use crate::scalar::{Algebra, Conjugate, Elementary, Field, Real, RealElementary};

impl<T: Algebra, const N: usize> Vector<T, N> {
    /// Summe aller Komponenten.
    #[inline]
    pub fn sum(self) -> T {
        self.0.into_iter().fold(T::ZERO, |acc, a| acc + a)
    }

    /// Komponentenweises Produkt.
    #[inline]
    pub fn component_mul(self, other: Self) -> Self {
        self.zip_map(other, |a, b| a * b)
    }

    /// Bilineares Skalarprodukt `Σ aᵢ·bᵢ` (ohne Konjugation).
    ///
    /// Für komplexe Vektoren ist meist [`inner`](Self::inner) gemeint.
    #[inline]
    pub fn dot(self, other: Self) -> T {
        self.component_mul(other).sum()
    }
}

impl<T: Algebra> Vector<T, 3> {
    /// Kreuzprodukt.
    #[inline]
    pub fn cross(self, other: Self) -> Self {
        let [a0, a1, a2] = self.0;
        let [b0, b1, b2] = other.0;
        Self([a1 * b2 - a2 * b1, a2 * b0 - a0 * b2, a0 * b1 - a1 * b0])
    }
}

impl<T: Conjugate, const N: usize> Vector<T, N> {
    /// Hermitesches Skalarprodukt `Σ conj(aᵢ)·bᵢ`.
    ///
    /// Für reelle Vektoren identisch mit [`dot`](Self::dot).
    #[inline]
    pub fn inner(self, other: Self) -> T {
        self.map(Conjugate::conj).dot(other)
    }

    /// Quadrat der euklidischen Norm.
    #[inline]
    pub fn norm_sqr(self) -> T::Real {
        self.0
            .into_iter()
            .fold(<T::Real as Algebra>::ZERO, |acc, a| acc + a.norm_sqr())
    }
}

impl<T: Conjugate, const N: usize> Vector<T, N>
where
    T::Real: RealElementary,
{
    /// Euklidische Norm (Länge).
    #[inline]
    pub fn norm(self) -> T::Real {
        self.norm_sqr().sqrt()
    }

    /// Auf Länge eins gebracht.
    ///
    /// Gibt `None` zurück, wenn die Länge null oder nicht endlich ist.
    pub fn normalize(self) -> Option<Self> {
        let norm = self.norm();
        if norm > <T::Real as Algebra>::ZERO && norm.is_finite() {
            let factor = norm.recip();
            Some(self.map(|a| a.scale(factor)))
        } else {
            None
        }
    }

    /// Euklidischer Abstand zweier Punkte.
    #[inline]
    pub fn distance(self, other: Self) -> T::Real {
        (self - other).norm()
    }
}
