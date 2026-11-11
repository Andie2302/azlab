//! Auswertung von Polynomen und Tschebyschow-Reihen aus Slices.

use crate::scalar::Algebra;

/// Wertet `c[0] + c[1]·x + c[2]·x² + …` mit dem Horner-Schema aus.
///
/// Für ein leeres Slice ist das Ergebnis null.
#[inline]
pub fn horner<T: Algebra>(coefficients: &[T], x: T) -> T {
    coefficients.iter().rev().fold(T::ZERO, |acc, &c| acc * x + c)
}

/// Wie [`horner`], liefert zusätzlich die Ableitung: `(Wert, Ableitung)`.
#[inline]
pub fn horner_with_derivative<T: Algebra>(coefficients: &[T], x: T) -> (T, T) {
    let mut value = T::ZERO;
    let mut derivative = T::ZERO;
    for &c in coefficients.iter().rev() {
        derivative = derivative * x + value;
        value = value * x + c;
    }
    (value, derivative)
}

/// Wertet die Tschebyschow-Reihe `c[0]·T₀(x) + c[1]·T₁(x) + …` mit dem
/// Clenshaw-Algorithmus aus.
///
/// Für ein leeres Slice ist das Ergebnis null.
pub fn clenshaw<T: Algebra>(coefficients: &[T], x: T) -> T {
    let Some((&first, rest)) = coefficients.split_first() else {
        return T::ZERO;
    };
    let two_x = x + x;
    let mut b1 = T::ZERO;
    let mut b2 = T::ZERO;
    for &c in rest.iter().rev() {
        let b0 = c + two_x * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    first + x * b1 - b2
}
