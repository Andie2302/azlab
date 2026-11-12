//! Sichten auf Matrizen, deren Daten in einem Puffer des Aufrufers liegen.
//!
//! Die Größe steht erst zur Laufzeit fest; es wird kein Heap gebraucht.
//! Gespeichert wird zeilenweise: Element `(i, j)` liegt an der Stelle
//! `i · cols + j`.

use super::LinalgError;
use super::decomp::{cholesky_in_place, left_divide, lu_in_place};
use crate::scalar::{Conjugate, Field, RealElementary};
use core::ops::{Index, IndexMut};

/// Position im Puffer; für ungültige Indizes eine, die sicher außerhalb liegt.
#[inline]
pub(super) fn flat_index(rows: usize, cols: usize, row: usize, col: usize) -> usize {
    if row < rows && col < cols { row * cols + col } else { usize::MAX }
}

/// Prüft, dass der Puffer genau `rows · cols` Elemente hat.
fn check_len(len: usize, rows: usize, cols: usize) -> Result<(), LinalgError> {
    match rows.checked_mul(cols) {
        Some(expected) if expected == len => Ok(()),
        _ => Err(LinalgError::DimensionMismatch),
    }
}

/// Lesende Sicht auf eine `rows × cols`-Matrix.
#[derive(Clone, Copy, Debug)]
pub struct MatrixView<'a, T> {
    data: &'a [T],
    rows: usize,
    cols: usize,
}

/// Schreibende Sicht auf eine `rows × cols`-Matrix.
#[derive(Debug)]
pub struct MatrixViewMut<'a, T> {
    data: &'a mut [T],
    rows: usize,
    cols: usize,
}

impl<'a, T> MatrixView<'a, T> {
    /// Sicht auf `data` als `rows × cols`-Matrix.
    ///
    /// Fehler, wenn `data` nicht genau `rows · cols` Elemente hat.
    pub fn new(data: &'a [T], rows: usize, cols: usize) -> Result<Self, LinalgError> {
        check_len(data.len(), rows, cols)?;
        Ok(Self { data, rows, cols })
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

    /// Die Elemente zeilenweise hintereinander.
    #[inline]
    pub fn data(&self) -> &'a [T] {
        self.data
    }

    /// Sicht ohne Längenprüfung; der Aufrufer garantiert
    /// `data.len() == rows · cols`.
    #[cfg(feature = "alloc")]
    #[inline]
    pub(super) fn from_parts(data: &'a [T], rows: usize, cols: usize) -> Self {
        Self { data, rows, cols }
    }
}

impl<'a, T> MatrixViewMut<'a, T> {
    /// Sicht auf `data` als `rows × cols`-Matrix.
    ///
    /// Fehler, wenn `data` nicht genau `rows · cols` Elemente hat.
    pub fn new(data: &'a mut [T], rows: usize, cols: usize) -> Result<Self, LinalgError> {
        check_len(data.len(), rows, cols)?;
        Ok(Self { data, rows, cols })
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

    /// Geprüfter schreibender Zugriff.
    #[inline]
    pub fn get_mut(&mut self, row: usize, col: usize) -> Option<&mut T> {
        self.data.get_mut(flat_index(self.rows, self.cols, row, col))
    }

    /// Sicht ohne Längenprüfung; der Aufrufer garantiert
    /// `data.len() == rows · cols`.
    #[cfg(feature = "alloc")]
    #[inline]
    pub(super) fn from_parts(data: &'a mut [T], rows: usize, cols: usize) -> Self {
        Self { data, rows, cols }
    }

    /// Lesende Sicht auf dieselben Daten.
    #[inline]
    pub fn as_view(&self) -> MatrixView<'_, T> {
        MatrixView { data: self.data, rows: self.rows, cols: self.cols }
    }
}

/// Zugriff `m[(zeile, spalte)]` wie bei Arrays: bricht bei falschem Index ab
/// (Panic). Für geprüften Zugriff gibt es [`MatrixView::get`].
impl<T> Index<(usize, usize)> for MatrixView<'_, T> {
    type Output = T;

    #[inline]
    fn index(&self, (row, col): (usize, usize)) -> &T {
        &self.data[flat_index(self.rows, self.cols, row, col)]
    }
}

/// Zugriff `m[(zeile, spalte)]` wie bei Arrays: bricht bei falschem Index ab
/// (Panic). Für geprüften Zugriff gibt es [`MatrixViewMut::get`].
impl<T> Index<(usize, usize)> for MatrixViewMut<'_, T> {
    type Output = T;

    #[inline]
    fn index(&self, (row, col): (usize, usize)) -> &T {
        &self.data[flat_index(self.rows, self.cols, row, col)]
    }
}

/// Schreibender Zugriff wie bei Arrays: bricht bei falschem Index ab (Panic).
/// Für geprüften Zugriff gibt es [`MatrixViewMut::get_mut`].
impl<T> IndexMut<(usize, usize)> for MatrixViewMut<'_, T> {
    #[inline]
    fn index_mut(&mut self, (row, col): (usize, usize)) -> &mut T {
        &mut self.data[flat_index(self.rows, self.cols, row, col)]
    }
}

impl<T: Field + Conjugate> MatrixViewMut<'_, T> {
    /// LU-Zerlegung mit Spaltenpivotisierung an Ort und Stelle.
    ///
    /// Danach steht `U` im oberen Dreieck und `L` ohne seine Einser-Diagonale
    /// darunter; `perm` (Länge `n`) nimmt die Zeilenpermutation auf.
    /// Rückgabe: `true` bei ungerader Zahl von Zeilentauschen. Zum Lösen
    /// dient [`MatrixView::lu_solve`].
    pub fn lu_in_place(&mut self, perm: &mut [usize]) -> Result<bool, LinalgError> {
        if self.rows != self.cols {
            return Err(LinalgError::DimensionMismatch);
        }
        let size = self.rows;
        lu_in_place(self, size, perm)
    }
}

impl<T: Field + Conjugate> MatrixViewMut<'_, T>
where
    T::Real: RealElementary,
{
    /// Cholesky-Zerlegung an Ort und Stelle; danach steht `L` im unteren
    /// Dreieck, das obere ist null.
    pub fn cholesky_in_place(&mut self) -> Result<(), LinalgError> {
        if self.rows != self.cols {
            return Err(LinalgError::DimensionMismatch);
        }
        let size = self.rows;
        cholesky_in_place(self, size)
    }
}

impl<T: Field + Conjugate> MatrixView<'_, T> {
    /// Löst `A·x = b` mit einer Matrix, die zuvor mit
    /// [`MatrixViewMut::lu_in_place`] zerlegt wurde; `perm` ist die dabei
    /// entstandene Permutation.
    // Die Schleifen brauchen den Index für Matrix und Vektor zugleich.
    #[allow(clippy::needless_range_loop)]
    pub fn lu_solve(&self, perm: &[usize], b: &[T], x: &mut [T]) -> Result<(), LinalgError> {
        let size = self.rows;
        if self.cols != size || perm.len() != size || b.len() != size || x.len() != size {
            return Err(LinalgError::DimensionMismatch);
        }
        // Vorwärtseinsetzen mit L (Einser-Diagonale).
        for i in 0..size {
            let Some(&start) = perm.get(i).and_then(|&source| b.get(source)) else {
                return Err(LinalgError::DimensionMismatch);
            };
            let mut sum = start;
            for j in 0..i {
                sum -= self[(i, j)] * x[j];
            }
            x[i] = sum;
        }
        // Rückwärtseinsetzen mit U.
        for i in (0..size).rev() {
            let mut sum = x[i];
            for j in i + 1..size {
                sum -= self[(i, j)] * x[j];
            }
            x[i] = left_divide(self[(i, i)], sum);
        }
        Ok(())
    }
}
