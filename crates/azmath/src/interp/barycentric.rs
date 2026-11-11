//! Polynominterpolation in baryzentrischer Form.

use super::error::InterpError;
use crate::scalar::Real;

/// Das eindeutige Polynom vom Grad höchstens `N − 1` durch `N` Punkte.
///
/// Die baryzentrische Form ist numerisch stabil und wertet in `O(N)` aus.
/// Für viele gleichabständige Punkte schwingt ein Interpolationspolynom
/// stark; dann sind Splines oder Tschebyschow-Stützstellen
/// ([`chebyshev_nodes`](super::chebyshev_nodes)) die bessere Wahl.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Barycentric<T, const N: usize> {
    xs: [T; N],
    ys: [T; N],
    weights: [T; N],
}

impl<T: Real, const N: usize> Barycentric<T, N> {
    /// Baut das Polynom durch die Punkte `(xs[i], ys[i])`.
    ///
    /// Fehler, wenn es keinen Punkt gibt oder zwei Stützstellen gleich sind.
    pub fn new(xs: [T; N], ys: [T; N]) -> Result<Self, InterpError> {
        if N == 0 {
            return Err(InterpError::TooFewPoints);
        }
        let mut weights = [T::ONE; N];
        for (j, weight) in weights.iter_mut().enumerate() {
            let mut product = T::ONE;
            for (k, &other) in xs.iter().enumerate() {
                if k != j {
                    product *= xs[j] - other;
                }
            }
            if product == T::ZERO || product.is_nan() {
                return Err(InterpError::NotIncreasing);
            }
            *weight = product.recip();
        }
        Ok(Self { xs, ys, weights })
    }

    /// Wert des Polynoms an der Stelle `x` (auch außerhalb der Stützstellen).
    pub fn eval(&self, x: T) -> T {
        let mut numerator = T::ZERO;
        let mut denominator = T::ZERO;
        for ((&node, &value), &weight) in self.xs.iter().zip(&self.ys).zip(&self.weights) {
            let distance = x - node;
            if distance == T::ZERO {
                return value;
            }
            let term = weight / distance;
            numerator += term * value;
            denominator += term;
        }
        numerator / denominator
    }
}
