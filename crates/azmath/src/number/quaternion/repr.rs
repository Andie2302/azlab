//! Der Typ [`Quaternion`].

use crate::scalar::Real;

/// Quaternion `w + x·i + y·j + z·k`.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Quaternion<T> {
    /// Skalarteil.
    pub w: T,
    /// Faktor vor `i`.
    pub x: T,
    /// Faktor vor `j`.
    pub y: T,
    /// Faktor vor `k`.
    pub z: T,
}

impl<T> Quaternion<T> {
    /// Erzeugt `w + x·i + y·j + z·k`.
    #[inline]
    pub const fn new(w: T, x: T, y: T, z: T) -> Self {
        Self { w, x, y, z }
    }
}

impl<T: Real> Quaternion<T> {
    /// Die Eins; als Rotation die Drehung um den Winkel null.
    pub const IDENTITY: Self = Self::new(T::ONE, T::ZERO, T::ZERO, T::ZERO);
    /// Imaginäre Einheit `i`.
    pub const I: Self = Self::new(T::ZERO, T::ONE, T::ZERO, T::ZERO);
    /// Imaginäre Einheit `j`.
    pub const J: Self = Self::new(T::ZERO, T::ZERO, T::ONE, T::ZERO);
    /// Imaginäre Einheit `k`.
    pub const K: Self = Self::new(T::ZERO, T::ZERO, T::ZERO, T::ONE);

    /// Reelle Zahl als Quaternion.
    #[inline]
    pub fn from_real(w: T) -> Self {
        Self::new(w, T::ZERO, T::ZERO, T::ZERO)
    }

    /// Quaternion aus Skalarteil und Vektorteil.
    #[inline]
    pub fn from_parts(w: T, [x, y, z]: [T; 3]) -> Self {
        Self::new(w, x, y, z)
    }

    /// Skalarteil `w`.
    #[inline]
    pub fn scalar(self) -> T {
        self.w
    }

    /// Vektorteil `[x, y, z]`.
    #[inline]
    pub fn vector(self) -> [T; 3] {
        [self.x, self.y, self.z]
    }

    /// Skalarprodukt der vier Komponenten.
    #[inline]
    pub fn dot(self, other: Self) -> T {
        self.w * other.w + self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Inverses Element `conj / norm_sqr`; für null ergeben sich NaN-Anteile.
    #[inline]
    pub fn inv(self) -> Self {
        let scale = self.dot(self).recip();
        Self::new(self.w * scale, -self.x * scale, -self.y * scale, -self.z * scale)
    }
}
