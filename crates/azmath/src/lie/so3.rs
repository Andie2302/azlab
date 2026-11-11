//! Die Drehgruppe des Raums, SO(3).

use super::tangent::So3Tangent;
use super::traits::LieGroup;
use crate::linalg::{Matrix, Vector};
use crate::number::Quaternion;
use crate::scalar::{Conjugate, RealElementary};

/// Drehung im Raum, gespeichert als Einheitsquaternion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct So3<T>(Quaternion<T>);

impl<T: RealElementary> So3<T> {
    /// Drehung aus einer Quaternion; sie wird normiert.
    ///
    /// Gibt `None` zurück, wenn ihre Norm null oder nicht endlich ist.
    #[inline]
    pub fn from_quaternion(quaternion: Quaternion<T>) -> Option<Self> {
        quaternion.normalize().map(Self)
    }

    /// Drehung um `axis` mit dem Winkel `angle` (Bogenmaß).
    ///
    /// Gibt `None` zurück, wenn die Achse die Länge null hat.
    #[inline]
    pub fn from_axis_angle(axis: Vector<T, 3>, angle: T) -> Option<Self> {
        Quaternion::from_axis_angle(axis.0, angle).map(Self)
    }

    /// Die Einheitsquaternion der Drehung.
    #[inline]
    pub fn quaternion(self) -> Quaternion<T> {
        self.0
    }

    /// Dreht den Vektor `v`.
    #[inline]
    pub fn rotate(self, v: Vector<T, 3>) -> Vector<T, 3> {
        Vector(self.0.rotate(v.0))
    }

    /// Die Drehmatrix; ihre Spalten sind die gedrehten Koordinatenachsen.
    pub fn matrix(self) -> Matrix<T, 3, 3> {
        let columns = [
            self.0.rotate([T::ONE, T::ZERO, T::ZERO]),
            self.0.rotate([T::ZERO, T::ONE, T::ZERO]),
            self.0.rotate([T::ZERO, T::ZERO, T::ONE]),
        ];
        Matrix::from_fn(|row, column| columns[column][row])
    }
}

impl<T: RealElementary> LieGroup for So3<T> {
    type Algebra = So3Tangent<T>;

    #[inline]
    fn identity() -> Self {
        Self(Quaternion::IDENTITY)
    }

    #[inline]
    fn compose(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }

    #[inline]
    fn inverse(self) -> Self {
        Self(self.0.conj())
    }

    /// Drehung um die Richtung von `phi` mit dem Winkel `|phi|`.
    fn exp(tangent: So3Tangent<T>) -> Self {
        let half = tangent.phi * T::from_f64(0.5);
        Self(Quaternion::from_parts(T::ZERO, half.0).exp())
    }

    /// Drehvektor mit Winkel in [0, π].
    fn log(self) -> So3Tangent<T> {
        // q und −q beschreiben dieselbe Drehung; die Form mit w >= 0 hat
        // den kleineren Winkel.
        let canonical = if self.0.w < T::ZERO { -self.0 } else { self.0 };
        let half = Vector(canonical.ln().vector());
        So3Tangent { phi: half * T::from_f64(2.0) }
    }

    #[inline]
    fn adjoint(self, tangent: So3Tangent<T>) -> So3Tangent<T> {
        So3Tangent { phi: self.rotate(tangent.phi) }
    }
}
