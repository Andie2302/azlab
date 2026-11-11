//! Lie-Gruppen und Lie-Algebren für Drehungen und starre Bewegungen.
//!
//! | Gruppe | Bedeutung | Algebra |
//! |--------|-----------|---------|
//! | [`So2`] | Drehungen in der Ebene | [`So2Tangent`] |
//! | [`So3`] | Drehungen im Raum | [`So3Tangent`] |
//! | [`Se2`] | starre Bewegungen in der Ebene | [`Se2Tangent`] |
//! | [`Se3`] | starre Bewegungen im Raum | [`Se3Tangent`] |
//!
//! Tangentialvektoren von [`Se2`] und [`Se3`] führen erst die Translation
//! (`rho`), dann die Drehung.

mod se2;
mod se3;
mod so2;
mod so3;
mod tangent;
mod traits;

#[cfg(test)]
mod tests;

pub use se2::Se2;
pub use se3::Se3;
pub use so2::So2;
pub use so3::So3;
pub use tangent::{Se2Tangent, Se3Tangent, So2Tangent, So3Tangent};
pub use traits::{LieAlgebra, LieGroup};
