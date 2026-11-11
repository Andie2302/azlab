//! Norm, Exponentialfunktion und Rotationen für [`Quaternion`].

use super::repr::Quaternion;
use crate::linalg::scaled_norm;
use crate::scalar::{Algebra, RealElementary};

fn cross<T: RealElementary>(a: [T; 3], b: [T; 3]) -> [T; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn length<T: RealElementary>([x, y, z]: [T; 3]) -> T {
    (x * x + y * y + z * z).sqrt()
}

impl<T: RealElementary> Quaternion<T> {
    /// Norm (Betrag), gegen Über- und Unterlauf geschützt.
    #[inline]
    pub fn norm(self) -> T {
        scaled_norm([self.w, self.x, self.y, self.z].into_iter())
    }

    /// Auf Norm eins gebracht.
    ///
    /// Gibt `None` zurück, wenn die Norm null oder nicht endlich ist.
    pub fn normalize(self) -> Option<Self> {
        let norm = self.norm();
        if norm > T::ZERO && norm.is_finite() {
            Some(self / norm)
        } else {
            None
        }
    }

    /// Exponentialfunktion.
    pub fn exp(self) -> Self {
        let angle = length(self.vector());
        let radius = self.w.exp();
        if angle == T::ZERO {
            return Self::from_real(radius);
        }
        let (sin, cos) = angle.sin_cos();
        let scale = radius * sin / angle;
        Self::new(radius * cos, self.x * scale, self.y * scale, self.z * scale)
    }

    /// Natürlicher Logarithmus (Hauptzweig, Winkel in [0, π]).
    ///
    /// Für negative reelle Zahlen ist die Drehachse nicht bestimmt; dann
    /// wird die `i`-Achse gewählt, also `ln(−1) = π·i`.
    pub fn ln(self) -> Self {
        let vector_norm = length(self.vector());
        let ln_norm = self.norm().ln();
        if vector_norm == T::ZERO {
            return if self.w < T::ZERO {
                Self::new(ln_norm, T::PI, T::ZERO, T::ZERO)
            } else {
                Self::from_real(ln_norm)
            };
        }
        let scale = vector_norm.atan2(self.w) / vector_norm;
        Self::new(ln_norm, self.x * scale, self.y * scale, self.z * scale)
    }

    /// Potenz mit reellem Exponenten, `exp(t · ln(self))`.
    ///
    /// Für Einheitsquaternionen ist das der Bruchteil `t` der Drehung.
    pub fn powf(self, t: T) -> Self {
        (self.ln() * t).exp()
    }

    /// Drehung um `axis` mit dem Winkel `angle` (Bogenmaß, Rechte-Hand-Regel).
    ///
    /// Die Achse muss nicht normiert sein. Gibt `None` zurück, wenn ihre
    /// Länge null oder nicht endlich ist.
    pub fn from_axis_angle(axis: [T; 3], angle: T) -> Option<Self> {
        let len = length(axis);
        if !(len > T::ZERO && len.is_finite()) {
            return None;
        }
        let (sin, cos) = (angle * T::from_f64(0.5)).sin_cos();
        let scale = sin / len;
        Some(Self::new(cos, axis[0] * scale, axis[1] * scale, axis[2] * scale))
    }

    /// Zerlegung einer Einheitsquaternion in `(Achse, Winkel)`.
    ///
    /// Der Winkel liegt in [0, 2π]. Für die Drehung um null wird die
    /// x-Achse zurückgegeben.
    pub fn to_axis_angle(self) -> ([T; 3], T) {
        let vector_norm = length(self.vector());
        if vector_norm == T::ZERO {
            return ([T::ONE, T::ZERO, T::ZERO], T::ZERO);
        }
        let half = vector_norm.atan2(self.w);
        (
            [self.x / vector_norm, self.y / vector_norm, self.z / vector_norm],
            half + half,
        )
    }

    /// Dreht den Vektor `v`. Setzt eine Einheitsquaternion voraus.
    pub fn rotate(self, v: [T; 3]) -> [T; 3] {
        // v' = v + w·t + u × t  mit  t = 2·(u × v)
        let u = self.vector();
        let t = cross(u, v).map(|c| c + c);
        let ut = cross(u, t);
        [
            v[0] + self.w * t[0] + ut[0],
            v[1] + self.w * t[1] + ut[1],
            v[2] + self.w * t[2] + ut[2],
        ]
    }

    /// Normierte lineare Interpolation auf dem kürzeren Weg.
    ///
    /// Schneller als [`slerp`](Self::slerp), aber mit nicht ganz konstanter
    /// Drehgeschwindigkeit. Heben sich beide Enden auf, wird `self`
    /// zurückgegeben.
    pub fn nlerp(self, other: Self, t: T) -> Self {
        let other = if self.dot(other) < T::ZERO { -other } else { other };
        (self * (T::ONE - t) + other * t).normalize().unwrap_or(self)
    }

    /// Sphärische Interpolation zwischen zwei Einheitsquaternionen auf dem
    /// kürzeren Weg, mit konstanter Drehgeschwindigkeit.
    pub fn slerp(self, other: Self, t: T) -> Self {
        let mut cos = self.dot(other);
        let mut other = other;
        if cos < T::ZERO {
            cos = -cos;
            other = -other;
        }
        // Fast gleiche Richtungen: sin(θ) wäre nahe null.
        if cos > T::from_f64(0.9995) {
            return self.nlerp(other, t);
        }
        let angle = cos.acos();
        let sin = angle.sin();
        let a = ((T::ONE - t) * angle).sin() / sin;
        let b = (t * angle).sin() / sin;
        self * a + other * b
    }
}

impl<T: RealElementary> Quaternion<T> {
    /// Quadratwurzel, `self.powf(0.5)`.
    #[inline]
    pub fn sqrt(self) -> Self {
        if self == Self::ZERO {
            return self;
        }
        self.powf(T::from_f64(0.5))
    }
}
