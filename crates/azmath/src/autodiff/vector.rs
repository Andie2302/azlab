//! Ableitungen von Funktionen mehrerer Veränderlicher.

use crate::linalg::{Matrix, Vector};
use crate::number::DualN;
use crate::scalar::Algebra;

/// Macht aus dem Punkt `x` die `N` unabhängigen Variablen: Komponente `i`
/// trägt den Wert `x[i]` und die Ableitung eins in Richtung `i`.
pub fn variables<T: Algebra, const N: usize>(x: Vector<T, N>) -> Vector<DualN<T, N>, N> {
    let mut index = 0;
    x.map(|value| {
        let mut eps = [T::ZERO; N];
        if let Some(slot) = eps.get_mut(index) {
            *slot = T::ONE;
        }
        index += 1;
        DualN::new(value, eps)
    })
}

/// Wert und Gradient von `f` an der Stelle `x`.
pub fn gradient<T, F, const N: usize>(f: F, x: Vector<T, N>) -> (T, Vector<T, N>)
where
    T: Algebra,
    F: FnOnce(Vector<DualN<T, N>, N>) -> DualN<T, N>,
{
    let value = f(variables(x));
    (value.re, Vector(value.eps))
}

/// Wert und Jacobi-Matrix von `f` an der Stelle `x`.
///
/// Zeile `i` der Matrix ist der Gradient der `i`-ten Komponente von `f`.
pub fn jacobian<T, F, const N: usize, const M: usize>(
    f: F,
    x: Vector<T, N>,
) -> (Vector<T, M>, Matrix<T, M, N>)
where
    T: Algebra,
    F: FnOnce(Vector<DualN<T, N>, N>) -> Vector<DualN<T, N>, M>,
{
    let value = f(variables(x));
    (value.map(|component| component.re), Matrix(value.0.map(|component| component.eps)))
}

/// Wert, Gradient und Hesse-Matrix von `f` an der Stelle `x`.
///
/// `f` wird dafür mit geschachtelten Dualzahlen aufgerufen.
pub fn hessian<T, F, const N: usize>(f: F, x: Vector<T, N>) -> (T, Vector<T, N>, Matrix<T, N, N>)
where
    T: Algebra,
    F: FnOnce(Vector<DualN<DualN<T, N>, N>, N>) -> DualN<DualN<T, N>, N>,
{
    let value = f(variables(variables(x)));
    (
        value.re.re,
        Vector(value.re.eps),
        Matrix(value.eps.map(|row| row.eps)),
    )
}
