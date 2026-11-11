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

/// Die Abbildung `x ↦ identity·x + scale·A·x` zu einem quadratischen
/// Operator `A`.
///
/// Damit entstehen aus einem Operator ohne Kopie verwandte Systeme, etwa
/// `−Δ` (für die Poisson-Gleichung) oder `I − dt·Δ` (für einen impliziten
/// Zeitschritt).
#[derive(Debug)]
pub struct AffineOperator<'a, A, T> {
    operator: &'a A,
    identity: T,
    scale: T,
}

impl<'a, A, T> AffineOperator<'a, A, T> {
    /// Erzeugt `identity·I + scale·operator`.
    #[inline]
    pub fn new(operator: &'a A, identity: T, scale: T) -> Self {
        Self { operator, identity, scale }
    }
}

impl<A: LinearOperator<T>, T: Algebra> LinearOperator<T> for AffineOperator<'_, A, T> {
    #[inline]
    fn rows(&self) -> usize {
        self.operator.rows()
    }

    #[inline]
    fn cols(&self) -> usize {
        self.operator.cols()
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if self.operator.rows() != self.operator.cols() {
            return Err(LinalgError::DimensionMismatch);
        }
        self.operator.apply(x, y)?;
        for (out, &input) in y.iter_mut().zip(x) {
            *out = self.identity * input + self.scale * *out;
        }
        Ok(())
    }
}
