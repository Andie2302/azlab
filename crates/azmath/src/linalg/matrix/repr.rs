//! Der Typ [`Matrix`] mit Zugriff und Umwandlungen.

use crate::linalg::Vector;
use crate::scalar::Algebra;
use core::ops::{Index, IndexMut};

/// Matrix mit `R` Zeilen und `C` Spalten auf dem Stack, zeilenweise
/// gespeichert: `m.0[zeile][spalte]`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matrix<T, const R: usize, const C: usize>(pub [[T; C]; R]);

impl<T, const R: usize, const C: usize> Matrix<T, R, C> {
    /// Erzeugt eine Matrix aus ihren Zeilen.
    #[inline]
    pub const fn new(rows: [[T; C]; R]) -> Self {
        Self(rows)
    }

    /// Erzeugt eine Matrix aus einer Funktion `(zeile, spalte) -> Wert`.
    #[inline]
    pub fn from_fn(mut f: impl FnMut(usize, usize) -> T) -> Self {
        Self(core::array::from_fn(|i| core::array::from_fn(|j| f(i, j))))
    }

    /// Anzahl der Zeilen.
    #[inline]
    pub const fn rows(&self) -> usize {
        R
    }

    /// Anzahl der Spalten.
    #[inline]
    pub const fn cols(&self) -> usize {
        C
    }

    /// Die Zeilen als verschachteltes Array.
    #[inline]
    pub const fn as_rows(&self) -> &[[T; C]; R] {
        &self.0
    }

    /// Geprüfter Zugriff: `None`, wenn Zeile oder Spalte außerhalb liegt.
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> Option<&T> {
        self.0.get(row)?.get(col)
    }

    /// Geprüfter schreibender Zugriff.
    #[inline]
    pub fn get_mut(&mut self, row: usize, col: usize) -> Option<&mut T> {
        self.0.get_mut(row)?.get_mut(col)
    }

    /// Wendet `f` auf jedes Element an.
    #[inline]
    pub fn map<U>(self, mut f: impl FnMut(T) -> U) -> Matrix<U, R, C> {
        Matrix(self.0.map(|row| row.map(&mut f)))
    }
}

impl<T: Copy, const R: usize, const C: usize> Matrix<T, R, C> {
    /// Zeile `index` als Vektor; `None` bei `index >= R`.
    #[inline]
    pub fn row(&self, index: usize) -> Option<Vector<T, C>> {
        self.0.get(index).map(|row| Vector(*row))
    }

    /// Spalte `index` als Vektor; `None` bei `index >= C`.
    #[inline]
    pub fn col(&self, index: usize) -> Option<Vector<T, R>> {
        if index < C {
            Some(Vector(self.0.map(|row| row[index])))
        } else {
            None
        }
    }

    /// Transponierte Matrix.
    #[inline]
    pub fn transpose(self) -> Matrix<T, C, R> {
        Matrix::from_fn(|i, j| self.0[j][i])
    }

    /// Verknüpft zwei Matrizen elementweise.
    #[inline]
    pub fn zip_map(self, other: Self, f: impl Fn(T, T) -> T) -> Self {
        Self::from_fn(|i, j| f(self.0[i][j], other.0[i][j]))
    }
}

impl<T: Algebra, const R: usize, const C: usize> Matrix<T, R, C> {
    /// Nullmatrix.
    pub const ZERO: Self = Self([[T::ZERO; C]; R]);
}

impl<T: Default + Copy, const R: usize, const C: usize> Default for Matrix<T, R, C> {
    #[inline]
    fn default() -> Self {
        Self([[T::default(); C]; R])
    }
}

impl<T, const R: usize, const C: usize> From<[[T; C]; R]> for Matrix<T, R, C> {
    #[inline]
    fn from(rows: [[T; C]; R]) -> Self {
        Self(rows)
    }
}

impl<T, const R: usize, const C: usize> From<Matrix<T, R, C>> for [[T; C]; R] {
    #[inline]
    fn from(m: Matrix<T, R, C>) -> Self {
        m.0
    }
}

/// Zugriff `m[(zeile, spalte)]` wie bei Arrays: bricht bei falschem Index ab
/// (Panic). Für geprüften Zugriff gibt es [`Matrix::get`].
impl<T, const R: usize, const C: usize> Index<(usize, usize)> for Matrix<T, R, C> {
    type Output = T;

    #[inline]
    fn index(&self, (row, col): (usize, usize)) -> &T {
        &self.0[row][col]
    }
}

/// Zugriff `m[(zeile, spalte)]` wie bei Arrays: bricht bei falschem Index ab
/// (Panic). Für geprüften Zugriff gibt es [`Matrix::get_mut`].
impl<T, const R: usize, const C: usize> IndexMut<(usize, usize)> for Matrix<T, R, C> {
    #[inline]
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut T {
        &mut self.0[row][col]
    }
}
