//! Die Traits [`Algebra`] und [`Commutative`].

use core::fmt::Debug;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Rechenoperationen `+ - *` mit Null und Eins.
///
/// Sagt bewusst weder Assoziativität noch Kommutativität der Multiplikation
/// zu, damit auch Oktonionen hineinpassen. Diese Zusagen kommen über
/// [`Ring`](super::Ring) und [`Commutative`].
pub trait Algebra:
    Copy
    + Debug
    + PartialEq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
{
    /// Neutrales Element der Addition.
    const ZERO: Self;
    /// Neutrales Element der Multiplikation.
    const ONE: Self;
}

/// Marker: Die Multiplikation ist kommutativ.
pub trait Commutative: Algebra {}
