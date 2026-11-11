//! Interpolation: aus Stützwerten eine Funktion machen.
//!
//! | Verfahren | Wofür |
//! |-----------|-------|
//! | [`lerp`], [`linear`] | stückweise linear, robust, stetig |
//! | [`CubicSpline`], [`natural_spline`], [`spline_eval`] | glatte Kurve durch viele Punkte |
//! | [`Barycentric`] | ein Polynom durch wenige Punkte |
//! | [`ChebyshevSeries`] | eine Funktion auf einem Intervall sehr genau annähern |

mod barycentric;
mod chebyshev;
mod error;
mod linear;
mod spline;

#[cfg(test)]
mod tests;

pub use barycentric::Barycentric;
pub use chebyshev::{ChebyshevSeries, chebyshev_nodes};
pub use error::InterpError;
pub use linear::{lerp, linear};
pub use spline::{CubicSpline, natural_spline, spline_derivative, spline_eval};
