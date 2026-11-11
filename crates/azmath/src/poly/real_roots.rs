//! Der Typ [`RealRoots`].

use crate::scalar::Real;

/// Bis zu `N` verschiedene reelle Nullstellen, aufsteigend sortiert, ohne
/// Heap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RealRoots<T, const N: usize> {
    roots: [T; N],
    len: usize,
}

impl<T: Real, const N: usize> RealRoots<T, N> {
    /// Leere Nullstellenmenge.
    #[inline]
    pub(super) fn empty() -> Self {
        Self { roots: [T::ZERO; N], len: 0 }
    }

    /// Nullstellenmenge aus einem Slice; Duplikate und Überzählige entfallen.
    pub(super) fn from_slice(values: &[T]) -> Self {
        let mut out = Self::empty();
        for &value in values {
            out.insert(value);
        }
        out
    }

    /// Fügt `value` sortiert ein. Exakt gleiche Werte, NaN und Werte über
    /// der Kapazität werden verworfen.
    pub(super) fn insert(&mut self, value: T) {
        if value.is_nan() || self.as_slice().contains(&value) {
            return;
        }
        let Some(slot) = self.roots.get_mut(self.len) else {
            return;
        };
        *slot = value;
        self.len += 1;
        // Einsortieren: nach vorn tauschen, solange der Vorgänger größer ist.
        let mut index = self.len - 1;
        while index > 0 && self.roots[index - 1] > self.roots[index] {
            self.roots.swap(index - 1, index);
            index -= 1;
        }
    }

    /// Die Nullstellen, aufsteigend.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        self.roots.get(..self.len).unwrap_or(&[])
    }

    /// Anzahl der Nullstellen.
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    /// `true`, wenn es keine reelle Nullstelle gibt.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Iterator über die Nullstellen, aufsteigend.
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
        self.as_slice().iter()
    }
}
