//! Zerlegungen: LU, Cholesky, QR, Eigenwerte, Singulärwerte.
//!
//! Die eigentlichen Algorithmen sind gegen `IndexMut<(usize, usize)>`
//! geschrieben und kennen den Typ [`Matrix`](crate::linalg::Matrix) nicht.
//! Andere Speicherformen (Sichten, dynamische Matrizen) können sie deshalb
//! später unverändert nutzen.
//!
//! Elementtypen: alles mit `Field + Conjugate`, also reelle und komplexe
//! Zahlen sowie Dualzahlen (Ableitungen laufen durch das Lösen hindurch).

mod cholesky;
mod eigen;
mod lu;
mod qr;
mod svd;
mod util;

#[cfg(test)]
mod spectral_tests;

#[cfg(test)]
mod tests;

pub(crate) use cholesky::cholesky_in_place;
pub(crate) use lu::lu_in_place;

pub use cholesky::Cholesky;
pub use eigen::SymmetricEigen;
pub use svd::Svd;
pub use lu::Lu;
pub use qr::Qr;
