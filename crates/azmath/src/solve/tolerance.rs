//! Abbruchkriterien.

use crate::scalar::Real;

/// Abbruchkriterien eines Lösers.
///
/// Ein Verfahren gilt als konvergiert, wenn der letzte Schritt höchstens
/// `abs_x + rel_x·|x|` groß war oder der Restfehler höchstens `abs_f` ist.
/// Der Standard zielt auf volle Maschinengenauigkeit.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct Tolerance<R> {
    /// Absolute Toleranz für die Schrittweite.
    pub abs_x: R,
    /// Relative Toleranz für die Schrittweite.
    pub rel_x: R,
    /// Absolute Toleranz für den Restfehler.
    pub abs_f: R,
    /// Höchstzahl der Iterationen.
    pub max_iterations: u32,
}

impl<R: Real> Default for Tolerance<R> {
    fn default() -> Self {
        Self {
            abs_x: R::EPSILON * R::EPSILON,
            rel_x: R::from_f64(4.0) * R::EPSILON,
            abs_f: R::ZERO,
            max_iterations: 200,
        }
    }
}

impl<R: Real> Tolerance<R> {
    /// Setzt die absolute Toleranz für die Schrittweite.
    #[inline]
    pub fn with_abs_x(mut self, abs_x: R) -> Self {
        self.abs_x = abs_x;
        self
    }

    /// Setzt die relative Toleranz für die Schrittweite.
    #[inline]
    pub fn with_rel_x(mut self, rel_x: R) -> Self {
        self.rel_x = rel_x;
        self
    }

    /// Setzt die absolute Toleranz für den Restfehler.
    #[inline]
    pub fn with_abs_f(mut self, abs_f: R) -> Self {
        self.abs_f = abs_f;
        self
    }

    /// Setzt die Höchstzahl der Iterationen.
    #[inline]
    pub fn with_max_iterations(mut self, max_iterations: u32) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// `true`, wenn ein Schritt der Größe `step` an der Stelle `x` klein
    /// genug ist.
    #[inline]
    pub(crate) fn x_converged(&self, step: R, x: R) -> bool {
        step.abs() <= self.abs_x + self.rel_x * x.abs()
    }

    /// `true`, wenn der Restfehler `residual` klein genug ist.
    #[inline]
    pub(crate) fn f_converged(&self, residual: R) -> bool {
        residual.abs() <= self.abs_f
    }
}
