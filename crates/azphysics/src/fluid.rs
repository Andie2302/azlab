//! Strömung: Navier-Stokes auf dem Gitter und Flüssigkeit aus Teilchen.
//!
//! - [`GridFluid`] löst die inkompressiblen Navier-Stokes-Gleichungen in der
//!   Ebene auf einem versetzten Gitter (*stable fluids*, J. Stam 1999):
//!   Transport entlang der Strömung, Zähigkeit, dann die Druckkorrektur, die
//!   das Geschwindigkeitsfeld quellenfrei macht. Die Druckgleichung löst das
//!   vorkonditionierte Verfahren der konjugierten Gradienten aus azmath. Geeignet für Rauch
//!   und für Umströmungen mit Einlass und Auslass.
//! - [`ParticleFluid`] rechnet eine Flüssigkeit mit freier Oberfläche aus
//!   Teilchen (*position based fluids*, Macklin und Müller 2013): Jedes
//!   Teilchen wird so verschoben, dass die Dichte nirgends über der
//!   Ruhedichte liegt.
//!
//! Beide brauchen das Feature `alloc`.

mod advect;
mod grid;
mod particles;
mod pressure;

#[cfg(all(test, any(feature = "std", feature = "libm")))]
mod tests;

pub use grid::{Edge, Edges, GridFluid, GridParams, GridStats};
pub use particles::{ParticleFluid, ParticleParams};

/// Punkt oder Richtung in der Ebene.
pub type Vec2<T> = azmath::Vector<T, 2>;
