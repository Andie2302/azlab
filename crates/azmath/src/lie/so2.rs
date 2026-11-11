//! Die Drehgruppe der Ebene, SO(2).

use super::tangent::So2Tangent;
use super::traits::LieGroup;
use crate::linalg::{Matrix, Vector};
use crate::scalar::RealElementary;

/// Drehung in der Ebene, gespeichert als Kosinus und Sinus des Winkels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct So2<T> {
    cos: T,
    sin: T,
}

impl<T: RealElementary> So2<T> {
    /// Drehung um `angle` (Bogenmaß, gegen den Uhrzeigersinn).
    #[inline]
    pub fn from_angle(angle: T) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self { cos, sin }
    }

    /// Drehwinkel im Bereich (−π, π].
    #[inline]
    pub fn angle(self) -> T {
        self.sin.atan2(self.cos)
    }

    /// Dreht den Vektor `v`.
    #[inline]
    pub fn rotate(self, v: Vector<T, 2>) -> Vector<T, 2> {
        let [x, y] = v.0;
        Vector([self.cos * x - self.sin * y, self.sin * x + self.cos * y])
    }

    /// Die Drehmatrix.
    #[inline]
    pub fn matrix(self) -> Matrix<T, 2, 2> {
        Matrix([[self.cos, -self.sin], [self.sin, self.cos]])
    }
}

impl<T: RealElementary> LieGroup for So2<T> {
    type Algebra = So2Tangent<T>;

    #[inline]
    fn identity() -> Self {
        Self { cos: T::ONE, sin: T::ZERO }
    }

    #[inline]
    fn compose(self, other: Self) -> Self {
        Self {
            cos: self.cos * other.cos - self.sin * other.sin,
            sin: self.sin * other.cos + self.cos * other.sin,
        }
    }

    #[inline]
    fn inverse(self) -> Self {
        Self { cos: self.cos, sin: -self.sin }
    }

    #[inline]
    fn exp(tangent: So2Tangent<T>) -> Self {
        Self::from_angle(tangent.angle)
    }

    #[inline]
    fn log(self) -> So2Tangent<T> {
        So2Tangent { angle: self.angle() }
    }

    /// Drehungen in der Ebene vertauschen: Die Adjungierte ändert nichts.
    #[inline]
    fn adjoint(self, tangent: So2Tangent<T>) -> So2Tangent<T> {
        tangent
    }
}
