//! QR-Zerlegung `A = Q·R` mit Householder-Spiegelungen.

use super::util::is_positive;
use crate::linalg::{LinalgError, Matrix, Vector};
use crate::scalar::{Algebra, Conjugate, Elementary, Field, RealElementary};
use core::ops::IndexMut;

/// Zerlegt die `rows × cols`-Matrix `a` an Ort und Stelle zu `R`.
///
/// `q` muss beim Aufruf die `rows × rows`-Einheitsmatrix sein und enthält
/// danach `Q`.
pub(crate) fn qr_in_place<T, A, Q>(a: &mut A, q: &mut Q, rows: usize, cols: usize)
where
    T: Field + Conjugate,
    T::Real: RealElementary,
    A: IndexMut<(usize, usize), Output = T>,
    Q: IndexMut<(usize, usize), Output = T>,
{
    let zero = <T::Real as Algebra>::ZERO;
    let two = <T::Real as Algebra>::ONE + <T::Real as Algebra>::ONE;

    for k in 0..rows.min(cols) {
        let mut column_sqr = zero;
        for i in k..rows {
            column_sqr += a[(i, k)].norm_sqr();
        }
        // Spalte ist schon null: nichts zu spiegeln. Die Prüfung vor der
        // Wurzel hält auch Ableitungen (Dualzahlen) endlich.
        if !is_positive(column_sqr) {
            continue;
        }
        let column_norm = column_sqr.sqrt();

        // Ziel: Spalte wird zu alpha·e₁ mit alpha = −phase·‖x‖. Die Phase
        // des Diagonalelements vermeidet Auslöschung und macht vᴴ·x reell.
        let head = a[(k, k)];
        let head_sqr = head.norm_sqr();
        let phase = if head_sqr > zero {
            head.scale(head_sqr.sqrt().recip())
        } else {
            T::ONE
        };
        let alpha = -phase.scale(column_norm);

        // Spiegelungsvektor v = x − alpha·e₁, abgelegt in der Spalte selbst.
        a[(k, k)] = head - alpha;
        let mut v_sqr = zero;
        for i in k..rows {
            v_sqr += a[(i, k)].norm_sqr();
        }
        let factor = two / v_sqr;

        // H = I − 2·v·vᴴ/(vᴴ·v) von links auf die restlichen Spalten.
        for j in k + 1..cols {
            let mut projection = T::ZERO;
            for i in k..rows {
                projection += a[(i, k)].conj() * a[(i, j)];
            }
            let weight = projection.scale(factor);
            for i in k..rows {
                let reduced = a[(i, j)] - a[(i, k)] * weight;
                a[(i, j)] = reduced;
            }
        }
        // Q ← Q·H.
        for r in 0..rows {
            let mut projection = T::ZERO;
            for i in k..rows {
                projection += q[(r, i)] * a[(i, k)];
            }
            let weight = projection.scale(factor);
            for i in k..rows {
                let reduced = q[(r, i)] - weight * a[(i, k)].conj();
                q[(r, i)] = reduced;
            }
        }

        a[(k, k)] = alpha;
        for i in k + 1..rows {
            a[(i, k)] = T::ZERO;
        }
    }
}

/// Ergebnis der QR-Zerlegung.
#[derive(Clone, Copy, Debug)]
pub struct Qr<T, const R: usize, const C: usize> {
    q: Matrix<T, R, R>,
    r: Matrix<T, R, C>,
}

impl<T: Field + Conjugate, const R: usize, const C: usize> Qr<T, R, C> {
    /// Der unitäre (reell: orthogonale) Faktor `Q`.
    pub fn q(&self) -> Matrix<T, R, R> {
        self.q
    }

    /// Der Faktor `R`: obere Dreiecks- bzw. Trapezform.
    pub fn r(&self) -> Matrix<T, R, C> {
        self.r
    }

    /// Löst `A·x = b` im Sinn der kleinsten Quadrate (`‖A·x − b‖` minimal).
    ///
    /// Für quadratische, reguläre Matrizen ist das die exakte Lösung.
    /// Fehler: [`LinalgError::DimensionMismatch`] bei weniger Zeilen als
    /// Spalten, [`LinalgError::Singular`] bei Rangabfall (ein
    /// Diagonalelement von `R` ist exakt null).
    pub fn solve(&self, b: Vector<T, R>) -> Result<Vector<T, C>, LinalgError> {
        if R < C {
            return Err(LinalgError::DimensionMismatch);
        }
        let projected = self.q.adjoint() * b;
        let mut x = [T::ZERO; C];
        for (i, row) in self.r.0.iter().enumerate().take(C).rev() {
            let diagonal = row[i];
            if !is_positive(diagonal.norm_sqr()) {
                return Err(LinalgError::Singular);
            }
            let rest = row[i + 1..]
                .iter()
                .zip(&x[i + 1..])
                .fold(projected.0[i], |acc, (&u, &xj)| acc - u * xj);
            x[i] = diagonal.recip() * rest;
        }
        Ok(Vector(x))
    }
}

impl<T: Field + Conjugate, const R: usize, const C: usize> Matrix<T, R, C>
where
    T::Real: RealElementary,
{
    /// QR-Zerlegung mit Householder-Spiegelungen.
    pub fn qr(self) -> Qr<T, R, C> {
        let mut r = self;
        let mut q = Matrix::<T, R, R>::identity();
        qr_in_place(&mut r, &mut q, R, C);
        Qr { q, r }
    }
}
