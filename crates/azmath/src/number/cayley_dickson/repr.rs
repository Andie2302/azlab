//! Der Typ [`CayleyDickson`].

use crate::scalar::{Conjugate, Field, RealElementary};

/// Paar `(a, b)` über der Algebra `A`, gelesen als `a + b·ℓ` mit einer neuen
/// imaginären Einheit `ℓ`.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct CayleyDickson<A> {
    /// Erster Anteil.
    pub a: A,
    /// Zweiter Anteil (Faktor vor der neuen imaginären Einheit).
    pub b: A,
}

/// Oktonionen: dreifache Verdopplung der reellen Zahlen.
pub type Octonion<T> = CayleyDickson<CayleyDickson<CayleyDickson<T>>>;

impl<A> CayleyDickson<A> {
    /// Erzeugt das Paar `(a, b)`.
    #[inline]
    pub const fn new(a: A, b: A) -> Self {
        Self { a, b }
    }
}

impl<A: Conjugate> CayleyDickson<A> {
    /// Inverses Element `conj / norm_sqr`.
    ///
    /// Für null ergeben sich NaN-Anteile. Ab der vierten Stufe (Sedenionen)
    /// gibt es Nullteiler; dort ist das Ergebnis kein echtes Inverses mehr.
    #[inline]
    pub fn inv(self) -> Self {
        self.conj().scale(self.norm_sqr().recip())
    }
}

impl<A: Conjugate> CayleyDickson<A>
where
    A::Real: RealElementary,
{
    /// Norm (Betrag).
    #[inline]
    pub fn norm(self) -> A::Real {
        use crate::scalar::Elementary;
        self.norm_sqr().sqrt()
    }
}
