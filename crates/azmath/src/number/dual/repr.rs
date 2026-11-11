//! Die Typen [`Dual`] und [`DualN`].

use crate::scalar::Algebra;

/// Dualzahl `re + eps·ε` mit `ε² = 0`: Wert und eine Ableitung.
#[derive(Clone, Copy, Debug, Default)]
pub struct Dual<T> {
    /// Wert.
    pub re: T,
    /// Ableitung.
    pub eps: T,
}

/// Dualzahl mit `N` Ableitungsrichtungen: Wert und Gradient.
#[derive(Clone, Copy, Debug)]
pub struct DualN<T, const N: usize> {
    /// Wert.
    pub re: T,
    /// Ableitungen nach den `N` Variablen.
    pub eps: [T; N],
}

impl<T> Dual<T> {
    /// Erzeugt eine Dualzahl aus Wert und Ableitung.
    #[inline]
    pub const fn new(re: T, eps: T) -> Self {
        Self { re, eps }
    }
}

impl<T: Algebra> Dual<T> {
    /// Konstante: Ableitung null.
    #[inline]
    pub fn constant(re: T) -> Self {
        Self::new(re, T::ZERO)
    }

    /// Variable, nach der abgeleitet wird: Ableitung eins.
    #[inline]
    pub fn variable(re: T) -> Self {
        Self::new(re, T::ONE)
    }
}

impl<T, const N: usize> DualN<T, N> {
    /// Erzeugt eine Dualzahl aus Wert und Gradient.
    #[inline]
    pub const fn new(re: T, eps: [T; N]) -> Self {
        Self { re, eps }
    }
}

impl<T: Algebra, const N: usize> DualN<T, N> {
    /// Konstante: alle Ableitungen null.
    #[inline]
    pub fn constant(re: T) -> Self {
        Self::new(re, [T::ZERO; N])
    }

    /// Die Variable Nummer `index` (ab 0): ihre eigene Ableitung ist eins.
    ///
    /// Gibt `None` zurück, wenn `index >= N`.
    #[inline]
    pub fn variable(re: T, index: usize) -> Option<Self> {
        let mut eps = [T::ZERO; N];
        *eps.get_mut(index)? = T::ONE;
        Some(Self::new(re, eps))
    }
}

impl<T: Default + Copy, const N: usize> Default for DualN<T, N> {
    #[inline]
    fn default() -> Self {
        Self::new(T::default(), [T::default(); N])
    }
}
