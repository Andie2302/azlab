//! Der Trait [`Field`].

use super::ring::Ring;
use core::ops::{Div, DivAssign};

/// Ring mit Division (Körper bzw. Schiefkörper).
pub trait Field: Ring + Div<Output = Self> + DivAssign {
    /// Wandelt eine `f64`-Konstante in diesen Typ um (ggf. mit Rundung).
    fn from_f64(value: f64) -> Self;

    /// Kehrwert `1 / self`.
    #[inline]
    fn recip(self) -> Self {
        Self::ONE / self
    }

    /// Ganzzahlige Potenz durch wiederholtes Quadrieren.
    ///
    /// Kommt ohne Backend aus und funktioniert für jeden `Field`-Typ.
    fn powi(self, n: i32) -> Self {
        let mut base = if n < 0 { self.recip() } else { self };
        let mut exp = n.unsigned_abs();
        let mut acc = Self::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                acc *= base;
            }
            exp >>= 1;
            if exp > 0 {
                base *= base;
            }
        }
        acc
    }
}
