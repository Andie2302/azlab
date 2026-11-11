//! Numerische Integration (Quadratur) von Funktionen einer Veränderlichen.
//!
//! | Verfahren | Wofür |
//! |-----------|-------|
//! | [`trapezoid_rule`], [`simpson`] | einfache Summenformeln mit fester Zahl von Teilintervallen |
//! | [`GaussLegendre`] | glatte Integranden: sehr genau mit wenigen Auswertungen |
//! | [`adaptive_simpson`] | Integranden mit lokal schwierigem Verlauf; liefert eine Fehlerschätzung |

mod adaptive;
mod gauss;
mod newton_cotes;

#[cfg(test)]
mod tests;

pub use adaptive::{QuadError, QuadOptions, Quadrature, adaptive_simpson};
pub use gauss::GaussLegendre;
pub use newton_cotes::{simpson, trapezoid_rule};
