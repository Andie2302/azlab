//! Cholesky-Zerlegung `A = L·Lᴴ` für hermitesch positiv definite Matrizen.

use super::util::is_positive;
use crate::linalg::{LinalgError, Matrix, Vector};
use crate::scalar::{Conjugate, Elementary, Field, Real, RealElementary};
use core::ops::IndexMut;

/// Zerlegt die `n × n`-Matrix `a` an Ort und Stelle.
///
/// Gelesen wird nur das untere Dreieck mit Diagonale. Danach steht dort `L`,
/// das obere Dreieck ist null.
pub(crate) fn cholesky_in_place<T, M>(a: &mut M, n: usize) -> Result<(), LinalgError>
where
    T: Field + Conjugate,
    T::Real: RealElementary,
    M: IndexMut<(usize, usize), Output = T>,
{
    for j in 0..n {
        let mut diagonal = a[(j, j)].real_part();
        for k in 0..j {
            diagonal -= a[(j, k)].norm_sqr();
        }
        if !(is_positive(diagonal) && diagonal.is_finite()) {
            return Err(LinalgError::NotPositiveDefinite);
        }
        let root = diagonal.sqrt();
        let inverse_root = root.recip();
        a[(j, j)] = T::ONE.scale(root);

        for i in j + 1..n {
            let mut sum = a[(i, j)];
            for k in 0..j {
                sum -= a[(i, k)] * a[(j, k)].conj();
            }
            a[(i, j)] = sum.scale(inverse_root);
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            a[(i, j)] = T::ZERO;
        }
    }
    Ok(())
}

/// Ergebnis der Cholesky-Zerlegung.
#[derive(Clone, Copy, Debug)]
pub struct Cholesky<T, const N: usize> {
    l: Matrix<T, N, N>,
}

impl<T: Field + Conjugate, const N: usize> Cholesky<T, N> {
    /// Der Faktor `L`: untere Dreiecksmatrix mit positiver Diagonale.
    pub fn l(&self) -> Matrix<T, N, N> {
        self.l
    }

    /// Löst `A·x = b`.
    pub fn solve(&self, b: Vector<T, N>) -> Vector<T, N> {
        let rows = &self.l.0;
        let mut x = b.0;
        // Vorwärtseinsetzen mit L.
        for (i, row) in rows.iter().enumerate() {
            let rest = row[..i]
                .iter()
                .zip(&x[..i])
                .fold(x[i], |acc, (&l, &xj)| acc - l * xj);
            x[i] = row[i].recip() * rest;
        }
        // Rückwärtseinsetzen mit Lᴴ: Zeile i von Lᴴ ist die konjugierte
        // Spalte i von L.
        for i in (0..N).rev() {
            let rest = rows[i + 1..]
                .iter()
                .zip(&x[i + 1..])
                .fold(x[i], |acc, (row, &xj)| acc - row[i].conj() * xj);
            x[i] = rows[i][i].recip() * rest;
        }
        Vector(x)
    }

    /// Determinante von `A`: Quadrat des Produkts der Diagonale von `L`.
    pub fn det(&self) -> T {
        let product = self.l.diagonal().0.into_iter().fold(T::ONE, |acc, d| acc * d);
        product * product
    }

    /// Inverse von `A`.
    pub fn inverse(&self) -> Matrix<T, N, N> {
        let identity = Matrix::<T, N, N>::identity();
        let columns: [Vector<T, N>; N] =
            core::array::from_fn(|k| self.solve(Vector(identity.0.map(|row| row[k]))));
        Matrix::from_fn(|i, k| columns[k].0[i])
    }
}

impl<T: Field + Conjugate, const N: usize> Matrix<T, N, N>
where
    T::Real: RealElementary,
{
    /// Cholesky-Zerlegung einer hermitesch positiv definiten Matrix.
    ///
    /// Gelesen wird nur das untere Dreieck; ob die Matrix wirklich
    /// hermitesch ist, wird nicht geprüft. Fehler
    /// [`LinalgError::NotPositiveDefinite`], wenn sie nicht positiv definit
    /// ist.
    pub fn cholesky(self) -> Result<Cholesky<T, N>, LinalgError> {
        let mut l = self;
        cholesky_in_place(&mut l, N)?;
        Ok(Cholesky { l })
    }
}
