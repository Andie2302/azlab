//! Die Gruppe der starren Bewegungen der Ebene, SE(2).

use super::so2::So2;
use super::tangent::Se2Tangent;
use super::traits::LieGroup;
use crate::linalg::Vector;
use crate::scalar::RealElementary;

/// Starre Bewegung in der Ebene: erst drehen, dann verschieben.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Se2<T> {
    /// Drehanteil.
    pub rotation: So2<T>,
    /// Verschiebung.
    pub translation: Vector<T, 2>,
}

/// Koeffizienten `(sin θ / θ, (1 − cos θ) / θ)` mit ihren Grenzwerten für
/// kleine Winkel.
fn coefficients<T: RealElementary>(theta: T) -> (T, T) {
    if theta * theta < T::EPSILON {
        (T::ONE, theta * T::from_f64(0.5))
    } else {
        (theta.sin() / theta, (T::ONE - theta.cos()) / theta)
    }
}

impl<T: RealElementary> Se2<T> {
    /// Bewegung aus Drehung und Verschiebung.
    #[inline]
    pub fn new(rotation: So2<T>, translation: Vector<T, 2>) -> Self {
        Self { rotation, translation }
    }

    /// Wendet die Bewegung auf den Punkt `point` an.
    #[inline]
    pub fn act(self, point: Vector<T, 2>) -> Vector<T, 2> {
        self.rotation.rotate(point) + self.translation
    }
}

impl<T: RealElementary> LieGroup for Se2<T> {
    type Algebra = Se2Tangent<T>;

    #[inline]
    fn identity() -> Self {
        Self::new(So2::identity(), Vector::ZERO)
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

    fn exp(tangent: Se2Tangent<T>) -> Self {
        let (a, b) = coefficients(tangent.theta);
        let [x, y] = tangent.rho.0;
        Self::new(
            So2::from_angle(tangent.theta),
            Vector([a * x - b * y, b * x + a * y]),
        )
    }

    fn log(self) -> Se2Tangent<T> {
        let theta = self.rotation.angle();
        let (a, b) = coefficients(theta);
        let scale = a * a + b * b;
        let [x, y] = self.translation.0;
        Se2Tangent {
            rho: Vector([(a * x + b * y) / scale, (a * y - b * x) / scale]),
            theta,
        }
    }

    fn adjoint(self, tangent: Se2Tangent<T>) -> Se2Tangent<T> {
        let [x, y] = self.translation.0;
        Se2Tangent {
            rho: self.rotation.rotate(tangent.rho) + Vector([y, -x]) * tangent.theta,
            theta: tangent.theta,
        }
    }
}
