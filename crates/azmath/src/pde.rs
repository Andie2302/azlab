//! Bausteine für partielle Differentialgleichungen auf regelmäßigen Gittern.
//!
//! - [`UniformGrid`] beschreibt ein Gitter in `D` Raumrichtungen und rechnet
//!   zwischen Gitterindex, Position im Puffer und Koordinaten um.
//! - [`Laplacian`] und [`Derivative`] sind Finite-Differenzen-Operatoren mit
//!   wählbaren Randbedingungen ([`Boundary`]). Sie implementieren
//!   [`LinearOperator`](crate::linalg::LinearOperator) und passen damit zu
//!   den iterativen Lösern.
//! - Linienmethode: Die Ortsdiskretisierung macht aus der partiellen eine
//!   gewöhnliche Differentialgleichung. Explizit integriert man sie mit
//!   [`rk4_slice`](crate::ode::rk4_slice), implizit mit [`theta_step`].
//! - [`poisson_multigrid_1d`] und [`poisson_multigrid_2d`] lösen die
//!   Poisson-Gleichung mit einem Mehrgitterverfahren in optimaler Zeit.
//!
//! Felder sind Slices des Aufrufers; der letzte Gitterindex läuft am
//! schnellsten. Nichts braucht Heap-Speicher.

mod grid;
mod multigrid;
mod stencil;
mod time;

#[cfg(test)]
mod tests;

pub use grid::UniformGrid;
pub use multigrid::{poisson_multigrid_1d, poisson_multigrid_2d};
pub use stencil::{Boundary, Derivative, Laplacian, Scheme};
pub use time::{LinearSolver, theta_step};
