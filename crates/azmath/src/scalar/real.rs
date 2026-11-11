//! Der Trait [`Real`].

use super::field::Field;

/// Reelle Zahl: geordneter Körper mit Betrag und Konstanten.
///
/// Alles in diesem Trait ist ohne Rechen-Backend verfügbar.
pub trait Real: Field + PartialOrd {
    /// Kreiszahl π.
    const PI: Self;
    /// Eulersche Zahl e.
    const E: Self;
    /// Maschinengenauigkeit: Abstand zwischen 1 und der nächstgrößeren Zahl.
    const EPSILON: Self;
    /// Positiv unendlich.
    const INFINITY: Self;
    /// Keine Zahl (NaN).
    const NAN: Self;

    /// Wandelt in `f64` um (ggf. mit Rundung).
    fn to_f64(self) -> f64;
    /// Betrag.
    fn abs(self) -> Self;
    /// Kleinere der beiden Zahlen; NaN wird ignoriert, wenn die andere Zahl gültig ist.
    fn min(self, other: Self) -> Self;
    /// Größere der beiden Zahlen; NaN wird ignoriert, wenn die andere Zahl gültig ist.
    fn max(self, other: Self) -> Self;
    /// Betrag von `self` mit dem Vorzeichen von `sign`.
    fn copysign(self, sign: Self) -> Self;
    /// `true`, wenn die Zahl weder unendlich noch NaN ist.
    fn is_finite(self) -> bool;
    /// `true`, wenn die Zahl NaN ist.
    fn is_nan(self) -> bool;
}
