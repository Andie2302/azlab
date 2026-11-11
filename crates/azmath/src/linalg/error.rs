//! Fehler der linearen Algebra.

use core::fmt;

/// Fehler beim Zerlegen und Lösen.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum LinalgError {
    /// Die Matrix ist singulär (ein Pivot ist exakt null).
    Singular,
    /// Die Matrix ist nicht positiv definit.
    NotPositiveDefinite,
    /// Die Größen passen nicht zusammen.
    DimensionMismatch,
}

impl fmt::Display for LinalgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Singular => "Matrix ist singulär",
            Self::NotPositiveDefinite => "Matrix ist nicht positiv definit",
            Self::DimensionMismatch => "Größen passen nicht zusammen",
        })
    }
}

impl core::error::Error for LinalgError {}
