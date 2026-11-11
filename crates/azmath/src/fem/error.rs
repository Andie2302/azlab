//! Fehler der Finite-Elemente-Bausteine.

use core::fmt;

/// Fehler beim Aufbau von Elementen, Netzen und Operatoren.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FemError {
    /// Ein Element ist entartet (Länge, Fläche oder Volumen null oder nicht
    /// endlich, oder ein Viereck ist verdreht).
    DegenerateElement,
    /// Das Netz ist fehlerhaft: Eine Knotennummer liegt außerhalb, oder die
    /// Größe der Elementmatrix passt nicht zur Knotenzahl.
    InvalidMesh,
    /// Die Länge eines übergebenen Puffers passt nicht.
    LengthMismatch,
}

impl fmt::Display for FemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::DegenerateElement => "entartetes Element",
            Self::InvalidMesh => "fehlerhaftes Netz",
            Self::LengthMismatch => "Länge passt nicht",
        })
    }
}

impl core::error::Error for FemError {}
