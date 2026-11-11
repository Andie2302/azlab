//! Fehler der Löser.

use super::solution::Solution;
use core::fmt;

/// Fehler eines Lösers.
///
/// Wo es einen letzten Stand gibt, steckt er im Fehler: Der Aufrufer kann
/// damit neu starten oder die Näherung bewusst weiterverwenden.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum SolveError<X, R> {
    /// Nicht innerhalb der erlaubten Iterationen konvergiert.
    NotConverged(Solution<X, R>),
    /// Die Funktionswerte an den Intervallenden haben nicht verschiedene
    /// Vorzeichen.
    NoBracket,
    /// Ein Zwischenwert war unendlich oder NaN.
    NonFinite(Solution<X, R>),
    /// Die Ableitung bzw. die Jacobi-Matrix war singulär.
    Singular(Solution<X, R>),
    /// Ungültige Eingabe, z. B. ein nicht endliches Intervall.
    InvalidInput,
    /// Der Ausgabepuffer des Aufrufers ist zu klein.
    CapacityExceeded,
}

impl<X, R> fmt::Display for SolveError<X, R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::NotConverged(_) => "Verfahren ist nicht konvergiert",
            Self::NoBracket => "kein Vorzeichenwechsel im Intervall",
            Self::NonFinite(_) => "Zwischenwert ist unendlich oder NaN",
            Self::Singular(_) => "Ableitung oder Jacobi-Matrix ist singulär",
            Self::InvalidInput => "ungültige Eingabe",
            Self::CapacityExceeded => "Ausgabepuffer ist zu klein",
        })
    }
}

impl<X: fmt::Debug, R: fmt::Debug> core::error::Error for SolveError<X, R> {}

/// Ergebnis der Löser für eine Unbekannte.
pub type ScalarResult<T> = Result<Solution<T, T>, SolveError<T, T>>;
