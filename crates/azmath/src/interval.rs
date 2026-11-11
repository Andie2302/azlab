//! Intervallarithmetik mit Rundung nach außen.
//!
//! Ein [`Interval`] steht für alle reellen Zahlen zwischen seinen Grenzen.
//! Jede Rechenoperation liefert ein Intervall, das das exakte Ergebnis für
//! **alle** Punkte der Eingangsintervalle enthält – Rundungsfehler
//! eingeschlossen. Damit lassen sich Aussagen beweisen, etwa „in diesem
//! Bereich liegt sicher keine Nullstelle“.
//!
//! [`Interval`] implementiert die Skalar-Traits bis [`Field`](crate::Field),
//! damit generische Algorithmen (Polynomauswertung, Matrizen) unverändert
//! darauf laufen. Die Rechengesetze gelten dabei nur als Einschluss: zum
//! Beispiel ist `x − x` nicht null, sondern ein Intervall um null.

mod ops;
mod repr;
mod traits;

#[cfg(test)]
mod tests;

pub use repr::Interval;
