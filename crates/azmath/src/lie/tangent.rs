//! Die Lie-Algebren so(2), so(3), se(2) und se(3).

use super::traits::LieAlgebra;
use crate::linalg::{Matrix, Vector};
use crate::scalar::Real;
use core::ops::{Add, Neg, Sub};

/// Komponentenweise Addition, Subtraktion und Negation.
macro_rules! impl_tangent_ops {
    ($ty:ident { $($field:ident),+ }) => {
        impl<T: Real> Add for $ty<T> {
            type Output = Self;

            #[inline]
            fn add(self, rhs: Self) -> Self {
                Self { $($field: self.$field + rhs.$field),+ }
            }
        }

        impl<T: Real> Sub for $ty<T> {
            type Output = Self;

            #[inline]
            fn sub(self, rhs: Self) -> Self {
                Self { $($field: self.$field - rhs.$field),+ }
            }
        }

        impl<T: Real> Neg for $ty<T> {
            type Output = Self;

            #[inline]
            fn neg(self) -> Self {
                Self { $($field: -self.$field),+ }
            }
        }
    };
}

/// Element von so(2): ein Drehwinkel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct So2Tangent<T> {
    /// Drehwinkel im Bogenmaß.
    pub angle: T,
}

/// Element von so(3): ein Drehvektor (Achse mal Winkel).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct So3Tangent<T> {
    /// Drehvektor.
    pub phi: Vector<T, 3>,
}

/// Element von se(2): Translation und Drehwinkel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Se2Tangent<T> {
    /// Translationsanteil.
    pub rho: Vector<T, 2>,
    /// Drehwinkel im Bogenmaß.
    pub theta: T,
}

/// Element von se(3): Translation und Drehvektor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Se3Tangent<T> {
    /// Translationsanteil.
    pub rho: Vector<T, 3>,
    /// Drehvektor.
    pub phi: Vector<T, 3>,
}

impl_tangent_ops!(So2Tangent { angle });
impl_tangent_ops!(So3Tangent { phi });
impl_tangent_ops!(Se2Tangent { rho, theta });
impl_tangent_ops!(Se3Tangent { rho, phi });

/// Vektor um 90° gegen den Uhrzeigersinn gedreht.
#[inline]
fn perpendicular<T: Real>(v: Vector<T, 2>) -> Vector<T, 2> {
    let [x, y] = v.0;
    Vector([-y, x])
}

impl<T: Real> So3Tangent<T> {
    /// Die schiefsymmetrische Matrix `K` mit `K·v = phi × v`.
    pub fn hat(self) -> Matrix<T, 3, 3> {
        let [x, y, z] = self.phi.0;
        Matrix([[T::ZERO, -z, y], [z, T::ZERO, -x], [-y, x, T::ZERO]])
    }
}

impl<T: Real> LieAlgebra for So2Tangent<T> {
    type Scalar = T;

    const ZERO: Self = Self { angle: T::ZERO };

    #[inline]
    fn scale(self, factor: T) -> Self {
        Self { angle: self.angle * factor }
    }

    /// Drehungen in der Ebene vertauschen: Die Klammer ist null.
    #[inline]
    fn bracket(self, _other: Self) -> Self {
        Self::ZERO
    }
}

impl<T: Real> LieAlgebra for So3Tangent<T> {
    type Scalar = T;

    const ZERO: Self = Self { phi: Vector::ZERO };

    #[inline]
    fn scale(self, factor: T) -> Self {
        Self { phi: self.phi * factor }
    }

    /// Die Klammer ist das Kreuzprodukt der Drehvektoren.
    #[inline]
    fn bracket(self, other: Self) -> Self {
        Self { phi: self.phi.cross(other.phi) }
    }
}

impl<T: Real> LieAlgebra for Se2Tangent<T> {
    type Scalar = T;

    const ZERO: Self = Self { rho: Vector::ZERO, theta: T::ZERO };

    #[inline]
    fn scale(self, factor: T) -> Self {
        Self { rho: self.rho * factor, theta: self.theta * factor }
    }

    fn bracket(self, other: Self) -> Self {
        Self {
            rho: perpendicular(other.rho) * self.theta - perpendicular(self.rho) * other.theta,
            theta: T::ZERO,
        }
    }
}

impl<T: Real> LieAlgebra for Se3Tangent<T> {
    type Scalar = T;

    const ZERO: Self = Self { rho: Vector::ZERO, phi: Vector::ZERO };

    #[inline]
    fn scale(self, factor: T) -> Self {
        Self { rho: self.rho * factor, phi: self.phi * factor }
    }

    fn bracket(self, other: Self) -> Self {
        Self {
            rho: self.phi.cross(other.rho) + self.rho.cross(other.phi),
            phi: self.phi.cross(other.phi),
        }
    }
}
