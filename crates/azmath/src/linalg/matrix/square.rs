//! Quadratische Matrizen.

use super::repr::Matrix;
use crate::linalg::Vector;
use crate::scalar::{Algebra, Ring};

impl<T: Algebra, const N: usize> Matrix<T, N, N> {
    /// Einheitsmatrix.
    #[inline]
    pub fn identity() -> Self {
        Self::from_diagonal(Vector::splat(T::ONE))
    }

    /// Diagonalmatrix mit den gegebenen Diagonalelementen.
    #[inline]
    pub fn from_diagonal(diagonal: Vector<T, N>) -> Self {
        let mut out = Self::ZERO;
        for (row, (index, value)) in out.0.iter_mut().zip(diagonal.0.into_iter().enumerate()) {
            if let Some(slot) = row.get_mut(index) {
                *slot = value;
            }
        }
        out
    }

    /// Die Diagonale als Vektor.
    #[inline]
    pub fn diagonal(self) -> Vector<T, N> {
        Vector(core::array::from_fn(|i| self.0[i][i]))
    }

    /// Spur: Summe der Diagonalelemente.
    #[inline]
    pub fn trace(self) -> T {
        self.diagonal().sum()
    }
}

/// Quadratische Matrizen bilden selbst eine Algebra mit der Einheitsmatrix
/// als Eins und können so überall eingesetzt werden, wo ein [`Ring`]
/// verlangt ist.
impl<T: Ring, const N: usize> Algebra for Matrix<T, N, N> {
    const ZERO: Self = Self([[T::ZERO; N]; N]);
    const ONE: Self = {
        let mut rows = [[T::ZERO; N]; N];
        let mut i = 0;
        while i < N {
            rows[i][i] = T::ONE;
            i += 1;
        }
        Self(rows)
    };
}

impl<T: Ring, const N: usize> Ring for Matrix<T, N, N> {}
