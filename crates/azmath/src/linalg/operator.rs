//! Der Trait [`LinearOperator`].

use super::view::{MatrixView, MatrixViewMut};
use super::{LinalgError, Matrix};
use crate::scalar::Algebra;

/// Lineare Abbildung, von der nur die Wirkung auf einen Vektor bekannt sein
/// muss.
///
/// Iterative Löser arbeiten ausschließlich gegen diesen Trait. Dichte und
/// dünnbesetzte Matrizen implementieren ihn; matrixfreie Verfahren (z. B.
/// ein Differenzenstern auf einem Gitter) implementieren ihn selbst.
pub trait LinearOperator<T> {
    /// Länge des Ergebnisvektors.
    fn rows(&self) -> usize;

    /// Länge des Eingangsvektors.
    fn cols(&self) -> usize;

    /// Berechnet `y = A·x`.
    ///
    /// Fehler, wenn `x` nicht `cols()` oder `y` nicht `rows()` Elemente hat.
    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError>;
}

/// Gemeinsame Umsetzung für zeilenweise gespeicherte dichte Daten.
fn apply_dense<T: Algebra>(
    data: &[T],
    rows: usize,
    cols: usize,
    x: &[T],
    y: &mut [T],
) -> Result<(), LinalgError> {
    if x.len() != cols || y.len() != rows {
        return Err(LinalgError::DimensionMismatch);
    }
    if cols == 0 {
        y.fill(T::ZERO);
        return Ok(());
    }
    for (row, out) in data.chunks_exact(cols).zip(y.iter_mut()) {
        *out = row.iter().zip(x).fold(T::ZERO, |acc, (&a, &b)| acc + a * b);
    }
    Ok(())
}

impl<T: Algebra, const R: usize, const C: usize> LinearOperator<T> for Matrix<T, R, C> {
    #[inline]
    fn rows(&self) -> usize {
        R
    }

    #[inline]
    fn cols(&self) -> usize {
        C
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if x.len() != C || y.len() != R {
            return Err(LinalgError::DimensionMismatch);
        }
        for (row, out) in self.0.iter().zip(y.iter_mut()) {
            *out = row.iter().zip(x).fold(T::ZERO, |acc, (&a, &b)| acc + a * b);
        }
        Ok(())
    }
}

impl<T: Algebra> LinearOperator<T> for MatrixView<'_, T> {
    #[inline]
    fn rows(&self) -> usize {
        MatrixView::rows(self)
    }

    #[inline]
    fn cols(&self) -> usize {
        MatrixView::cols(self)
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        apply_dense(self.data(), MatrixView::rows(self), MatrixView::cols(self), x, y)
    }
}

impl<T: Algebra> LinearOperator<T> for MatrixViewMut<'_, T> {
    #[inline]
    fn rows(&self) -> usize {
        MatrixViewMut::rows(self)
    }

    #[inline]
    fn cols(&self) -> usize {
        MatrixViewMut::cols(self)
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        self.as_view().apply(x, y)
    }
}
