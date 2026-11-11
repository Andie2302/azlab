//! Die Gruppe der starren Bewegungen des Raums, SE(3).

use super::so3::So3;
use super::tangent::{Se3Tangent, So3Tangent};
use super::traits::LieGroup;
use crate::linalg::Vector;
use crate::scalar::RealElementary;

/// Starre Bewegung im Raum: erst drehen, dann verschieben.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Se3<T> {
    /// Drehanteil.
    pub rotation: So3<T>,
    /// Verschiebung.
    pub translation: Vector<T, 3>,
}

impl<T: RealElementary> Se3<T> {
    /// Bewegung aus Drehung und Verschiebung.
    #[inline]
    pub fn new(rotation: So3<T>, translation: Vector<T, 3>) -> Self {
        Self { rotation, translation }
    }

    /// Wendet die Bewegung auf den Punkt `point` an.
    #[inline]
    pub fn act(self, point: Vector<T, 3>) -> Vector<T, 3> {
        self.rotation.rotate(point) + self.translation
    }
}

impl<T: RealElementary> LieGroup for Se3<T> {
    type Algebra = Se3Tangent<T>;

    #[inline]
    fn identity() -> Self {
        Self::new(So3::identity(), Vector::ZERO)
    }

    #[inline]
    fn compose(self, other: Self) -> Self {
        Self::new(self.rotation.compose(other.rotation), self.act(other.translation))
    }

    #[inline]
    fn inverse(self) -> Self {
        let rotation = self.rotation.inverse();
        Self::new(rotation, -rotation.rotate(self.translation))
    }

    /// Translation `V·rho` mit
    /// `V = I + (1 − cos θ)/θ² · K + (θ − sin θ)/θ³ · K²`, `K = hat(phi)`.
    fn exp(tangent: Se3Tangent<T>) -> Self {
        let phi = tangent.phi;
        let theta_sqr = phi.dot(phi);
        let (a, b) = if theta_sqr < T::EPSILON {
            (T::from_f64(0.5), T::ONE / T::from_f64(6.0))
        } else {
            let theta = theta_sqr.sqrt();
            (
                (T::ONE - theta.cos()) / theta_sqr,
                (theta - theta.sin()) / (theta_sqr * theta),
            )
        };
        let k = phi.cross(tangent.rho);
        let kk = phi.cross(k);
        Self::new(So3::exp(So3Tangent { phi }), tangent.rho + k * a + kk * b)
    }

    /// Translation `V⁻¹·t` mit
    /// `V⁻¹ = I − K/2 + (1/θ² − (1 + cos θ)/(2θ·sin θ)) · K²`.
    ///
    /// Nahe einer halben Drehung (θ ≈ π) wird der letzte Koeffizient ungenau.
    fn log(self) -> Se3Tangent<T> {
        let phi = self.rotation.log().phi;
        let theta_sqr = phi.dot(phi);
        let c = if theta_sqr < T::EPSILON {
            T::ONE / T::from_f64(12.0)
        } else {
            let theta = theta_sqr.sqrt();
            T::ONE / theta_sqr
                - (T::ONE + theta.cos()) / (T::from_f64(2.0) * theta * theta.sin())
        };
        let k = phi.cross(self.translation);
        let kk = phi.cross(k);
        Se3Tangent {
            rho: self.translation - k * T::from_f64(0.5) + kk * c,
            phi,
        }
    }

    fn adjoint(self, tangent: Se3Tangent<T>) -> Se3Tangent<T> {
        let phi = self.rotation.rotate(tangent.phi);
        Se3Tangent {
            rho: self.rotation.rotate(tangent.rho) + self.translation.cross(phi),
            phi,
        }
    }
}
