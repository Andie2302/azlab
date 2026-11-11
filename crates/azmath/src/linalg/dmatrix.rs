//! Besitzende Matrix mit Größe zur Laufzeit (Feature `alloc`).

use super::view::{MatrixView, MatrixViewMut, flat_index};
use super::{LinalgError, LinearOperator};
use crate::scalar::{Algebra, Conjugate, Field};
use alloc::vec;
use alloc::vec::Vec;
use core::ops::{Index, IndexMut};

/// Matrix auf dem Heap, zeilenweise gespeichert; die Größe steht erst zur
/// Laufzeit fest.
#[derive(Clone, Debug, PartialEq)]
pub struct DMatrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}

impl<T: Algebra> DMatrix<T> {
    /// Nullmatrix. Fehler, wenn `rows · cols` nicht darstellbar ist.
    pub fn zeros(rows: usize, cols: usize) -> Result<Self, LinalgError> {
        let len = rows.checked_mul(cols).ok_or(LinalgError::DimensionMismatch)?;
        Ok(Self { data: vec![T::ZERO; len], rows, cols })
    }

    /// Einheitsmatrix der Größe `size × size`.
    pub fn identity(size: usize) -> Result<Self, LinalgError> {
        let mut out = Self::zeros(size, size)?;
        for index in 0..size {
            out[(index, index)] = T::ONE;
        }
        Ok(out)
    }

    /// Matrix aus zeilenweise abgelegten Elementen.
    ///
    /// Fehler, wenn `data` nicht genau `rows · cols` Elemente hat.
    pub fn from_vec(data: Vec<T>, rows: usize, cols: usize) -> Result<Self, LinalgError> {
        if rows.checked_mul(cols) == Some(data.len()) {
            Ok(Self { data, rows, cols })
        } else {
            Err(LinalgError::DimensionMismatch)
        }
    }

    /// Matrix aus einer Funktion `(zeile, spalte) -> Wert`.
    pub fn from_fn(
        rows: usize,
        cols: usize,
        mut f: impl FnMut(usize, usize) -> T,
    ) -> Result<Self, LinalgError> {
        let mut out = Self::zeros(rows, cols)?;
        for row in 0..rows {
            for col in 0..cols {
                out[(row, col)] = f(row, col);
            }
        }
        Ok(out)
    }

    /// Anzahl der Zeilen.
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Anzahl der Spalten.
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Geprüfter Zugriff: `None`, wenn Zeile oder Spalte außerhalb liegt.
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> Option<&T> {
        self.data.get(flat_index(self.rows, self.cols, row, col))
    }

    /// Lesende Sicht.
    #[inline]
    pub fn view(&self) -> MatrixView<'_, T> {
        MatrixView::from_parts(&self.data, self.rows, self.cols)
    }

    /// Schreibende Sicht.
    #[inline]
    pub fn view_mut(&mut self) -> MatrixViewMut<'_, T> {
        MatrixViewMut::from_parts(&mut self.data, self.rows, self.cols)
    }

    /// Transponierte Matrix.
    pub fn transpose(&self) -> Self {
        let mut out = Self { data: vec![T::ZERO; self.data.len()], rows: self.cols, cols: self.rows };
        for row in 0..self.rows {
            for col in 0..self.cols {
                out[(col, row)] = self[(row, col)];
            }
        }
        out
    }

    /// Matrixprodukt. Fehler, wenn die Größen nicht zusammenpassen.
    // Kein `Mul`-Operator, weil das Produkt scheitern kann.
    #[allow(clippy::should_implement_trait)]
    pub fn mul(&self, other: &Self) -> Result<Self, LinalgError> {
        if self.cols != other.rows {
            return Err(LinalgError::DimensionMismatch);
        }
        let mut out = Self::zeros(self.rows, other.cols)?;
        for row in 0..self.rows {
            for col in 0..other.cols {
                let mut sum = T::ZERO;
                for inner in 0..self.cols {
                    sum += self[(row, inner)] * other[(inner, col)];
                }
                out[(row, col)] = sum;
            }
        }
        Ok(out)
    }

    /// Matrix mal Vektor. Fehler, wenn `x` nicht `cols` Elemente hat.
    pub fn mul_vec(&self, x: &[T]) -> Result<Vec<T>, LinalgError> {
        let mut y = vec![T::ZERO; self.rows];
        self.view().apply(x, &mut y)?;
        Ok(y)
    }
}

impl<T: Field + Conjugate> DMatrix<T> {
    /// Löst `self·x = b` über die LU-Zerlegung einer Kopie.
    pub fn solve(&self, b: &[T]) -> Result<Vec<T>, LinalgError> {
        let mut factors = self.clone();
        let mut perm = vec![0; self.rows];
        factors.view_mut().lu_in_place(&mut perm)?;
        let mut x = vec![T::ZERO; self.rows];
        factors.view().lu_solve(&perm, b, &mut x)?;
        Ok(x)
    }
}

/// Zugriff `m[(zeile, spalte)]` wie bei Arrays: bricht bei falschem Index ab
/// (Panic). Für geprüften Zugriff gibt es [`DMatrix::get`].
impl<T> Index<(usize, usize)> for DMatrix<T> {
    type Output = T;

    #[inline]
    fn index(&self, (row, col): (usize, usize)) -> &T {
        &self.data[flat_index(self.rows, self.cols, row, col)]
    }
}

/// Schreibender Zugriff wie bei Arrays: bricht bei falschem Index ab (Panic).
impl<T> IndexMut<(usize, usize)> for DMatrix<T> {
    #[inline]
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut T {
        &mut self.data[flat_index(self.rows, self.cols, row, col)]
    }
}

impl<T: Algebra> LinearOperator<T> for DMatrix<T> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn cols(&self) -> usize {
        self.cols
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        self.view().apply(x, y)
    }
}
