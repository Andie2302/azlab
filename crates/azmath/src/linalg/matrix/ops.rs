//! Rechenoperatoren und Normen für [`Matrix`].

use super::repr::Matrix;
use crate::linalg::Vector;
use crate::scalar::{Algebra, Conjugate, Elementary, Field, RealElementary};
use core::ops::{Add, Div, DivAssign, Mul, MulAssign, Neg, Sub};

impl<T: Algebra, const R: usize, const C: usize> Add for Matrix<T, R, C> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self {
        self.zip_map(rhs, |a, b| a + b)
    }
}

impl<T: Algebra, const R: usize, const C: usize> Sub for Matrix<T, R, C> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        self.zip_map(rhs, |a, b| a - b)
    }
}

impl<T: Algebra, const R: usize, const C: usize> Neg for Matrix<T, R, C> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        self.map(|a| -a)
    }
}

/// Multiplikation jedes Elements mit einem Skalar (von rechts).
impl<T: Algebra, const R: usize, const C: usize> Mul<T> for Matrix<T, R, C> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self {
        self.map(|a| a * rhs)
    }
}

/// Division jedes Elements durch einen Skalar.
impl<T: Field, const R: usize, const C: usize> Div<T> for Matrix<T, R, C> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self {
        self.map(|a| a / rhs)
    }
}

/// Matrixprodukt; die Größen prüft der Compiler.
impl<T: Algebra, const R: usize, const C: usize, const K: usize> Mul<Matrix<T, C, K>>
    for Matrix<T, R, C>
{
    type Output = Matrix<T, R, K>;

    fn mul(self, rhs: Matrix<T, C, K>) -> Matrix<T, R, K> {
        Matrix::from_fn(|i, k| {
            self.0[i]
                .iter()
                .zip(&rhs.0)
                .fold(T::ZERO, |acc, (&a, row)| acc + a * row[k])
        })
    }
}

/// Matrix mal Spaltenvektor.
impl<T: Algebra, const R: usize, const C: usize> Mul<Vector<T, C>> for Matrix<T, R, C> {
    type Output = Vector<T, R>;

    fn mul(self, rhs: Vector<T, C>) -> Vector<T, R> {
        Vector(self.0.map(|row| Vector(row).dot(rhs)))
    }
}

impl_op_assign!([T: Algebra, const R: usize, const C: usize] Matrix<T, R, C>, AddAssign, add_assign, +);
impl_op_assign!([T: Algebra, const R: usize, const C: usize] Matrix<T, R, C>, SubAssign, sub_assign, -);
impl_op_assign!([T: Algebra, const N: usize] Matrix<T, N, N>, MulAssign, mul_assign, *);

impl<T: Algebra, const R: usize, const C: usize> MulAssign<T> for Matrix<T, R, C> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

impl<T: Field, const R: usize, const C: usize> DivAssign<T> for Matrix<T, R, C> {
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}

// `T * Matrix<T, R, C>` lässt sich wegen der Kohärenzregeln nicht generisch
// schreiben, deshalb je primitivem Typ.
macro_rules! impl_left_scalar {
    ($($t:ty),*) => {
        $(
            impl<const R: usize, const C: usize> Mul<Matrix<$t, R, C>> for $t {
                type Output = Matrix<$t, R, C>;

                #[inline]
                fn mul(self, rhs: Matrix<$t, R, C>) -> Matrix<$t, R, C> {
                    rhs * self
                }
            }
        )*
    };
}

impl_left_scalar!(f32, f64);

impl<T: Algebra, const N: usize> Vector<T, N> {
    /// Äußeres Produkt `self · otherᵀ` (ohne Konjugation): Element
    /// `(i, j)` ist `self[i] · other[j]`.
    #[inline]
    pub fn outer<const M: usize>(self, other: Vector<T, M>) -> Matrix<T, N, M> {
        Matrix(self.0.map(|a| other.0.map(|b| a * b)))
    }
}

impl<T: Conjugate, const R: usize, const C: usize> Matrix<T, R, C> {
    /// Adjungierte (konjugiert transponierte) Matrix.
    #[inline]
    pub fn adjoint(self) -> Matrix<T, C, R> {
        self.transpose().map(Conjugate::conj)
    }

    /// Quadrat der Frobenius-Norm: Summe der Betragsquadrate aller Elemente.
    #[inline]
    pub fn norm_sqr(self) -> T::Real {
        self.0
            .into_iter()
            .fold(<T::Real as Algebra>::ZERO, |acc, row| acc + Vector(row).norm_sqr())
    }
}

impl<T: Conjugate, const R: usize, const C: usize> Matrix<T, R, C>
where
    T::Real: RealElementary,
{
    /// Frobenius-Norm.
    #[inline]
    pub fn norm(self) -> T::Real {
        self.norm_sqr().sqrt()
    }
}
