//! Der Typ [`Vector`] mit Zugriff und Umwandlungen.

use crate::scalar::Algebra;
use core::ops::{Index, IndexMut};

/// Vektor mit `N` Komponenten auf dem Stack.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector<T, const N: usize>(pub [T; N]);

impl<T, const N: usize> Vector<T, N> {
    /// Erzeugt einen Vektor aus seinen Komponenten.
    #[inline]
    pub const fn new(data: [T; N]) -> Self {
        Self(data)
    }

    /// Die Komponenten als Array.
    #[inline]
    pub const fn as_array(&self) -> &[T; N] {
        &self.0
    }

    /// Die Komponenten als Slice.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }

    /// Iterator über die Komponenten.
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
        self.0.iter()
    }

    /// Geprüfter Zugriff: `None` bei `index >= N`.
    #[inline]
    pub fn get(&self, index: usize) -> Option<&T> {
        self.0.get(index)
    }

    /// Geprüfter schreibender Zugriff: `None` bei `index >= N`.
    #[inline]
    pub fn get_mut(&mut self, index: usize) -> Option<&mut T> {
        self.0.get_mut(index)
    }

    /// Wendet `f` auf jede Komponente an.
    #[inline]
    pub fn map<U>(self, f: impl FnMut(T) -> U) -> Vector<U, N> {
        Vector(self.0.map(f))
    }
}

impl<T: Copy, const N: usize> Vector<T, N> {
    /// Verknüpft zwei Vektoren komponentenweise.
    #[inline]
    pub fn zip_map(self, other: Self, f: impl Fn(T, T) -> T) -> Self {
        let mut out = self;
        for (slot, rhs) in out.0.iter_mut().zip(other.0) {
            *slot = f(*slot, rhs);
        }
        out
    }
}

impl<T: Algebra, const N: usize> Vector<T, N> {
    /// Nullvektor.
    pub const ZERO: Self = Self([T::ZERO; N]);

    /// Vektor, dessen Komponenten alle gleich `value` sind.
    #[inline]
    pub fn splat(value: T) -> Self {
        Self([value; N])
    }

    /// Einheitsvektor in Richtung `index`; `None` bei `index >= N`.
    #[inline]
    pub fn unit(index: usize) -> Option<Self> {
        let mut out = Self::ZERO;
        *out.0.get_mut(index)? = T::ONE;
        Some(out)
    }
}

impl<T: Default + Copy, const N: usize> Default for Vector<T, N> {
    #[inline]
    fn default() -> Self {
        Self([T::default(); N])
    }
}

impl<T, const N: usize> From<[T; N]> for Vector<T, N> {
    #[inline]
    fn from(data: [T; N]) -> Self {
        Self(data)
    }
}

impl<T, const N: usize> From<Vector<T, N>> for [T; N] {
    #[inline]
    fn from(v: Vector<T, N>) -> Self {
        v.0
    }
}

/// Zugriff wie bei Arrays: bricht bei `index >= N` ab (Panic).
/// Für geprüften Zugriff gibt es [`Vector::get`].
impl<T, const N: usize> Index<usize> for Vector<T, N> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &T {
        &self.0[index]
    }
}

/// Zugriff wie bei Arrays: bricht bei `index >= N` ab (Panic).
/// Für geprüften Zugriff gibt es [`Vector::get_mut`].
impl<T, const N: usize> IndexMut<usize> for Vector<T, N> {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut T {
        &mut self.0[index]
    }
}

impl<'a, T, const N: usize> IntoIterator for &'a Vector<T, N> {
    type Item = &'a T;
    type IntoIter = core::slice::Iter<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}
