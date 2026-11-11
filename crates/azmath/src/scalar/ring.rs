//! Der Trait [`Ring`].

use super::algebra::Algebra;

/// Marker: Die Multiplikation ist assoziativ (Ring mit Eins).
///
/// Kommutativ muss sie nicht sein (Quaternionen).
pub trait Ring: Algebra {}
