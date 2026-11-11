//! Zusammensetzen der Gesamtmatrix im CSR-Format (Feature `alloc`).

use super::error::FemError;
use crate::linalg::{CsrMatrix, LinalgError, Matrix};
use crate::scalar::Real;
use alloc::vec::Vec;

/// Dünnbesetzte Matrix im CSR-Format mit eigenem Speicher.
#[derive(Clone, Debug, PartialEq)]
pub struct CsrOwned<T> {
    size: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<T>,
}

impl<T> CsrOwned<T> {
    /// Sicht für Löser und Operatoren.
    pub fn view(&self) -> Result<CsrMatrix<'_, T>, LinalgError> {
        CsrMatrix::new(self.size, self.size, &self.row_ptr, &self.col_idx, &self.values)
    }

    /// Anzahl der gespeicherten Elemente.
    #[inline]
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
}

/// Setzt die Gesamtmatrix aus den Elementmatrizen zusammen.
///
/// Bedeutung von `D`, `V`, `M` und der Nummerierung der Freiheitsgrade wie
/// bei [`ElementOperator`](super::ElementOperator). Randbedingungen werden
/// hier nicht eingearbeitet.
pub fn assemble_csr<T, F, const D: usize, const V: usize, const M: usize>(
    nodes: &[[T; D]],
    elements: &[[usize; V]],
    element_matrix: F,
) -> Result<CsrOwned<T>, FemError>
where
    T: Real,
    F: Fn(&[[T; D]; V]) -> Result<Matrix<T, M, M>, FemError>,
{
    if V == 0 || M == 0 || M % V != 0 {
        return Err(FemError::InvalidMesh);
    }
    let dofs_per_node = M / V;
    let size = nodes.len().checked_mul(dofs_per_node).ok_or(FemError::InvalidMesh)?;

    // Alle Beiträge als (Zeile, Spalte, Wert) sammeln und sortieren.
    let mut entries: Vec<(usize, usize, T)> = Vec::new();
    for element in elements {
        let mut corners = [[T::ZERO; D]; V];
        for (corner, &index) in corners.iter_mut().zip(element) {
            *corner = *nodes.get(index).ok_or(FemError::InvalidMesh)?;
        }
        let matrix = element_matrix(&corners)?;
        let dof = |local: usize| element[local / dofs_per_node] * dofs_per_node + local % dofs_per_node;
        for row in 0..M {
            for col in 0..M {
                entries.push((dof(row), dof(col), matrix[(row, col)]));
            }
        }
    }
    entries.sort_unstable_by_key(|&(row, col, _)| (row, col));

    // Gleiche Stellen zusammenfassen und Zeilenanfänge festhalten.
    let mut row_ptr = Vec::with_capacity(size + 1);
    let mut col_idx = Vec::new();
    let mut values: Vec<T> = Vec::new();
    row_ptr.push(0);
    let mut previous: Option<(usize, usize)> = None;
    for (row, col, value) in entries {
        if previous == Some((row, col)) {
            if let Some(last) = values.last_mut() {
                *last += value;
            }
            continue;
        }
        while row_ptr.len() <= row {
            row_ptr.push(values.len());
        }
        col_idx.push(col);
        values.push(value);
        previous = Some((row, col));
    }
    while row_ptr.len() <= size {
        row_ptr.push(values.len());
    }
    Ok(CsrOwned { size, row_ptr, col_idx, values })
}
