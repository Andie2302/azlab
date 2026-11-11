//! Iterative Löser für große lineare Gleichungssysteme `A·x = b`.
//!
//! Die Löser kennen die Matrix nur als
//! [`LinearOperator`](crate::linalg::LinearOperator) und brauchen keinen
//! Heap: Arbeitsspeicher stellt der Aufrufer bereit.
//!
//! Ein Verfahren gilt als konvergiert, wenn
//! `‖b − A·x‖ <= abs_f + rel_x · ‖b‖` gilt. Die Standardwerte von
//! [`Tolerance`](crate::solve::Tolerance) zielen auf volle
//! Maschinengenauigkeit, die iterative Verfahren oft nicht erreichen; in der
//! Regel setzt man `rel_x` passend zum Problem.

mod bicgstab;
mod cg;
mod gmres;
mod util;

#[cfg(test)]
mod tests;

pub use bicgstab::bicgstab;
pub use cg::conjugate_gradient;
pub use gmres::{gmres, gmres_work_len};
pub use util::LinearResult;
