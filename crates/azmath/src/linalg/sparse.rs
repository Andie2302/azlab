//! Dünnbesetzte Matrizen im CSR-Format als Sicht auf Puffer des Aufrufers.

use super::operator::LinearOperator;
use super::LinalgError;
use crate::scalar::Algebra;

/// Dünnbesetzte Matrix im CSR-Format (*compressed sparse row*).
///
/// - `values`: die gespeicherten Elemente, zeilenweise hintereinander,
/// - `col_idx`: zu jedem Element seine Spalte,
/// - `row_ptr`: Zeile `i` belegt die Positionen `row_ptr[i]..row_ptr[i + 1]`.
#[derive(Clone, Copy, Debug)]
pub struct CsrMatrix<'a, T> {
    rows: usize,
    cols: usize,
    row_ptr: &'a [usize],
    col_idx: &'a [usize],
    values: &'a [T],
}

impl<'a, T> CsrMatrix<'a, T> {
    /// Erzeugt die Sicht und prüft die Struktur.
    ///
    /// Fehler, wenn `row_ptr` nicht `rows + 1` Einträge hat, nicht bei null
    /// beginnt, nicht monoton wächst oder nicht bei der Zahl der Elemente
    /// endet, wenn `col_idx` und `values` verschieden lang sind oder wenn
    /// eine Spaltennummer außerhalb liegt.
    pub fn new(
        rows: usize,
        cols: usize,
        row_ptr: &'a [usize],
        col_idx: &'a [usize],
        values: &'a [T],
    ) -> Result<Self, LinalgError> {
        let consistent = row_ptr.len() == rows.wrapping_add(1)
            && row_ptr.first() == Some(&0)
            && row_ptr.last() == Some(&values.len())
            && row_ptr.windows(2).all(|pair| pair[0] <= pair[1])
            && col_idx.len() == values.len()
            && col_idx.iter().all(|&col| col < cols);
        if consistent {
            Ok(Self { rows, cols, row_ptr, col_idx, values })
        } else {
            Err(LinalgError::DimensionMismatch)
        }
    }

    /// Anzahl der gespeicherten Elemente.
    #[inline]
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Spalten und Werte der Zeile `row`; leer, wenn sie außerhalb liegt.
    fn row(&self, row: usize) -> (&'a [usize], &'a [T]) {
        let range = match (self.row_ptr.get(row), self.row_ptr.get(row.wrapping_add(1))) {
            (Some(&start), Some(&end)) => start..end,
            _ => 0..0,
        };
        (
            self.col_idx.get(range.clone()).unwrap_or(&[]),
            self.values.get(range).unwrap_or(&[]),
        )
    }
}

impl<T: Algebra> CsrMatrix<'_, T> {
    /// Element `(row, col)`; null, wenn es nicht gespeichert ist, `None`,
    /// wenn die Stelle außerhalb liegt. Mehrfach gespeicherte Stellen werden
    /// addiert.
    pub fn get(&self, row: usize, col: usize) -> Option<T> {
        if row >= self.rows || col >= self.cols {
            return None;
        }
        let (columns, values) = self.row(row);
        Some(
            columns
                .iter()
                .zip(values)
                .filter(|&(&c, _)| c == col)
                .fold(T::ZERO, |acc, (_, &value)| acc + value),
        )
    }
}

impl<T: Algebra> LinearOperator<T> for CsrMatrix<'_, T> {
    #[inline]
    fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    fn cols(&self) -> usize {
        self.cols
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if x.len() != self.cols || y.len() != self.rows {
            return Err(LinalgError::DimensionMismatch);
        }
        for (row, out) in y.iter_mut().enumerate() {
            let (columns, values) = self.row(row);
            *out = columns.iter().zip(values).fold(T::ZERO, |acc, (&col, &value)| {
                // Die Spaltennummern wurden beim Erzeugen geprüft.
                acc + value * x.get(col).copied().unwrap_or(T::ZERO)
            });
        }
        Ok(())
    }
}
