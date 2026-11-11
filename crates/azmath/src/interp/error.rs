//! Fehler der Interpolation.

use core::fmt;

/// Fehler beim Aufbau einer Interpolation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum InterpError {
    /// Zu wenige Stützstellen.
    TooFewPoints,
    /// Die Stützstellen sind nicht streng aufsteigend bzw. nicht verschieden.
    NotIncreasing,
    /// Die Längen der übergebenen Puffer passen nicht zusammen.
    LengthMismatch,
}

impl fmt::Display for InterpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::TooFewPoints => "zu wenige Stützstellen",
            Self::NotIncreasing => "Stützstellen nicht streng aufsteigend",
            Self::LengthMismatch => "Längen passen nicht zusammen",
        })
    }
}

impl core::error::Error for InterpError {}
