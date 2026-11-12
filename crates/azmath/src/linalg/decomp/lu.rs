//! LU-Zerlegung mit Spaltenpivotisierung: `P·A = L·U`.

use crate::linalg::{LinalgError, Matrix, Vector};
use super::util::{is_positive, left_divide};
use crate::scalar::{Commutative, Conjugate, Field};
use core::ops::IndexMut;

/// Zerlegt die `n × n`-Matrix `a` an Ort und Stelle.
///
/// Danach steht `U` im oberen Dreieck (mit Diagonale) und `L` ohne seine
/// Einser-Diagonale darunter. `perm[i]` nennt die ursprüngliche Zeile, die
/// jetzt an Stelle `i` steht. Rückgabe: `true` bei ungerader Zahl von
/// Zeilentauschen.
pub(crate) fn lu_in_place<T, M>(a: &mut M, n: usize, perm: &mut [usize]) -> Result<bool, LinalgError>
where
    T: Field + Conjugate,
    M: IndexMut<(usize, usize), Output = T>,
{
    if perm.len() != n {
        return Err(LinalgError::DimensionMismatch);
    }
    for (index, slot) in perm.iter_mut().enumerate() {
        *slot = index;
    }

    let mut odd = false;
    for k in 0..n {
        // Pivot: betragsgrößtes Element der Spalte ab der Diagonale.
        let mut pivot = k;
        // Größe ohne Quadrieren: Sehr kleine oder sehr große Elemente
        // laufen so weder unter noch über.
        let mut best = a[(k, k)].max_abs();
        for i in k + 1..n {
            let magnitude = a[(i, k)].max_abs();
            if magnitude > best {
                best = magnitude;
                pivot = i;
            }
        }
        // Auch NaN gilt als nicht zerlegbar.
        if !is_positive(best) {
            return Err(LinalgError::Singular);
        }
        if pivot != k {
            for j in 0..n {
                let upper = a[(k, j)];
                a[(k, j)] = a[(pivot, j)];
                a[(pivot, j)] = upper;
            }
            perm.swap(k, pivot);
            odd = !odd;
        }

        let diagonal = a[(k, k)];
        for i in k + 1..n {
            let factor = a[(i, k)] / diagonal;
            a[(i, k)] = factor;
            for j in k + 1..n {
                let reduced = a[(i, j)] - factor * a[(k, j)];
                a[(i, j)] = reduced;
            }
        }
    }
    Ok(odd)
}

/// Ergebnis der LU-Zerlegung einer quadratischen Matrix.
///
/// Einmal zerlegen, dann beliebig viele rechte Seiten lösen.
#[derive(Clone, Copy, Debug)]
pub struct Lu<T, const N: usize> {
    lu: Matrix<T, N, N>,
    perm: [usize; N],
    odd: bool,
}

impl<T: Field + Conjugate, const N: usize> Lu<T, N> {
    /// Löst `A·x = b`.
    pub fn solve(&self, b: Vector<T, N>) -> Vector<T, N> {
        let rows = &self.lu.0;
        let mut x = [T::ZERO; N];
        // Vorwärtseinsetzen mit L (Einser-Diagonale).
        for (i, row) in rows.iter().enumerate() {
            let start = b.0[self.perm[i]];
            x[i] = row[..i]
                .iter()
                .zip(&x[..i])
                .fold(start, |acc, (&l, &xj)| acc - l * xj);
        }
        // Rückwärtseinsetzen mit U.
        for (i, row) in rows.iter().enumerate().rev() {
            let rest = row[i + 1..]
                .iter()
                .zip(&x[i + 1..])
                .fold(x[i], |acc, (&u, &xj)| acc - u * xj);
            x[i] = left_divide(row[i], rest);
        }
        Vector(x)
    }

    /// Löst `A·X = B` für mehrere rechte Seiten (die Spalten von `B`).
    pub fn solve_matrix<const K: usize>(&self, b: Matrix<T, N, K>) -> Matrix<T, N, K> {
        let columns: [Vector<T, N>; K] =
            core::array::from_fn(|k| self.solve(Vector(b.0.map(|row| row[k]))));
        Matrix::from_fn(|i, k| columns[k].0[i])
    }

    /// Inverse von `A`.
    pub fn inverse(&self) -> Matrix<T, N, N> {
        self.solve_matrix(Matrix::identity())
    }

    /// Der Faktor `L`: untere Dreiecksmatrix mit Einsen auf der Diagonale.
    pub fn l(&self) -> Matrix<T, N, N> {
        Matrix::from_fn(|i, j| match i.cmp(&j) {
            core::cmp::Ordering::Greater => self.lu.0[i][j],
            core::cmp::Ordering::Equal => T::ONE,
            core::cmp::Ordering::Less => T::ZERO,
        })
    }

    /// Der Faktor `U`: obere Dreiecksmatrix.
    pub fn u(&self) -> Matrix<T, N, N> {
        Matrix::from_fn(|i, j| if i <= j { self.lu.0[i][j] } else { T::ZERO })
    }

    /// Zeilenpermutation: Zeile `i` von `P·A` ist Zeile `permutation()[i]`
    /// von `A`.
    pub fn permutation(&self) -> [usize; N] {
        self.perm
    }
}

impl<T: Field + Conjugate, const N: usize> Matrix<T, N, N> {
    /// LU-Zerlegung mit Spaltenpivotisierung.
    ///
    /// Fehler [`LinalgError::Singular`], wenn ein Pivot exakt null ist.
    pub fn lu(self) -> Result<Lu<T, N>, LinalgError> {
        let mut lu = self;
        let mut perm = [0; N];
        let odd = lu_in_place(&mut lu, N, &mut perm)?;
        Ok(Lu { lu, perm, odd })
    }

    /// Inverse Matrix.
    pub fn inverse(self) -> Result<Self, LinalgError> {
        Ok(self.lu()?.inverse())
    }

    /// Löst `self·x = b`.
    ///
    /// Für mehrere rechte Seiten zur selben Matrix lohnt sich
    /// [`lu`](Self::lu) und dann [`Lu::solve`].
    pub fn solve(self, b: Vector<T, N>) -> Result<Vector<T, N>, LinalgError> {
        Ok(self.lu()?.solve(b))
    }
}

/// Die Determinante setzt eine kommutative Multiplikation voraus; für
/// Quaternionen etwa gibt es sie in dieser Form nicht.
impl<T: Field + Conjugate + Commutative, const N: usize> Lu<T, N> {
    /// Determinante von `A`.
    pub fn det(&self) -> T {
        let product = self
            .lu
            .diagonal()
            .0
            .into_iter()
            .fold(T::ONE, |acc, d| acc * d);
        if self.odd { -product } else { product }
    }
}

impl<T: Field + Conjugate + Commutative, const N: usize> Matrix<T, N, N> {
    /// Determinante. Für singuläre Matrizen null.
    ///
    /// Bis 3×3 über die geschlossene Formel, sonst über die LU-Zerlegung.
    pub fn det(self) -> T {
        let m = &self.0;
        match N {
            0 => T::ONE,
            1 => m[0][0],
            2 => m[0][0] * m[1][1] - m[0][1] * m[1][0],
            3 => {
                m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
            }
            _ => self.lu().map_or(T::ZERO, |lu| lu.det()),
        }
    }
}
