//! Zustand und Fehler der Differentialgleichungslöser.

use crate::linalg::Vector;
use core::fmt;

/// Zustand einer Integration: erreichte Zeit, Lösung dort und Schrittzähler.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OdeState<T, const N: usize> {
    /// Erreichte Zeit.
    pub t: T,
    /// Lösung zur Zeit `t`.
    pub y: Vector<T, N>,
    /// Anzahl der angenommenen Schritte.
    pub steps: u32,
    /// Anzahl der verworfenen Schritte (nur bei gesteuerter Schrittweite).
    pub rejected: u32,
}

/// Fehler eines Differentialgleichungslösers.
///
/// Wo es einen letzten Zustand gibt, steckt er im Fehler.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum OdeError<T, const N: usize> {
    /// Ungültige Eingabe, z. B. nicht endliche Zeiten oder null Schritte.
    InvalidInput,
    /// Die Lösung wurde unendlich oder NaN.
    NonFinite(OdeState<T, N>),
    /// Die erlaubte Zahl von Schritten reichte nicht.
    TooManySteps(OdeState<T, N>),
    /// Die Schrittweite wurde zu klein, um voranzukommen.
    StepTooSmall(OdeState<T, N>),
    /// Das Newton-Verfahren eines impliziten Schritts ist gescheitert.
    SolverFailed(OdeState<T, N>),
}

impl<T, const N: usize> fmt::Display for OdeError<T, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidInput => "ungültige Eingabe",
            Self::NonFinite(_) => "Lösung ist unendlich oder NaN",
            Self::TooManySteps(_) => "zu viele Schritte",
            Self::StepTooSmall(_) => "Schrittweite zu klein",
            Self::SolverFailed(_) => "impliziter Schritt nicht lösbar",
        })
    }
}

impl<T: fmt::Debug, const N: usize> core::error::Error for OdeError<T, N> {}

/// Ergebnis der Differentialgleichungslöser.
pub type OdeResult<T, const N: usize> = Result<OdeState<T, N>, OdeError<T, N>>;
