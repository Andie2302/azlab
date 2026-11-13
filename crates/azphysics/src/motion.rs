//! Fehler beim Wiederherstellen von Bewegungsdaten.
use core::fmt;

/// Ungültiger Bewegungszustand; bei einem Fehler bleibt die Simulation unverändert.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MotionError {
    /// Die Zahl der Punkte oder Geschwindigkeiten passt nicht zur Szene.
    DimensionMismatch,
    /// Eine Koordinate oder Geschwindigkeit ist nicht endlich.
    NonFinite,
    /// Ein Flüssigkeitspunkt liegt außerhalb des Beckens oder im Hindernis.
    InvalidPosition,
    /// Ein festgehaltener Punkt soll eine Geschwindigkeit erhalten.
    FixedParticle,
}

impl fmt::Display for MotionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::DimensionMismatch => "Anzahl der Bewegungsdaten passt nicht",
            Self::NonFinite => "nicht endliche Bewegungsdaten",
            Self::InvalidPosition => "Punkt außerhalb des gültigen Flüssigkeitsgebiets",
            Self::FixedParticle => "Geschwindigkeit für festgehaltenen Punkt",
        })
    }
}
impl core::error::Error for MotionError {}
