//! Der Trait [`Ring`].

use core::fmt::Debug;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Ring mit Eins: Addition, Subtraktion, Multiplikation.
///
/// Die Multiplikation muss assoziativ, aber nicht kommutativ sein
/// (Quaternionen).
pub trait Ring:
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
