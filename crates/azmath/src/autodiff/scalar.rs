//! Ableitungen von Funktionen einer Veränderlichen.

use crate::number::Dual;
use crate::scalar::Algebra;

/// Wert und erste Ableitung von `f` an der Stelle `x`: `(f, f')`.
pub fn derivative<T, F>(f: F, x: T) -> (T, T)
where
    T: Algebra,
    F: FnOnce(Dual<T>) -> Dual<T>,
{
    let value = f(Dual::variable(x));
    (value.re, value.eps)
}

/// Wert, erste und zweite Ableitung von `f` an der Stelle `x`:
/// `(f, f', f'')`.
pub fn second_derivative<T, F>(f: F, x: T) -> (T, T, T)
where
    T: Algebra,
    F: FnOnce(Dual<Dual<T>>) -> Dual<Dual<T>>,
{
    let value = f(Dual::variable(Dual::variable(x)));
    (value.re.re, value.re.eps, value.eps.eps)
}
