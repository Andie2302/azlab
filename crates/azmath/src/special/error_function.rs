//! Fehlerfunktion.

use super::gamma::{gamma_p, gamma_q};
use crate::scalar::RealElementary;

/// Fehlerfunktion `erf(x) = (2/√π)·∫₀ˣ e^(−t²) dt`.
///
/// Berechnet über `erf(x) = P(1/2, x²)` für `x >= 0`; ungerade Funktion.
/// Für große `|x|`, deren Quadrat nicht mehr darstellbar ist, gilt der
/// Grenzwert ±1.
pub fn erf<T: RealElementary>(x: T) -> T {
    let value = gamma_p(T::from_f64(0.5), x * x);
    if x < T::ZERO { -value } else { value }
}

/// Komplementäre Fehlerfunktion `erfc(x) = 1 − erf(x)`, ohne Auslöschung
/// für große `x`. Für große `|x|`, deren Quadrat nicht mehr darstellbar ist,
/// gelten die Grenzwerte 0 und 2.
pub fn erfc<T: RealElementary>(x: T) -> T {
    let half = T::from_f64(0.5);
    if x < T::ZERO {
        T::ONE + gamma_p(half, x * x)
    } else {
        gamma_q(half, x * x)
    }
}
