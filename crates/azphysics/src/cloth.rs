//! Gewebe: Tuch aus Teilchen und Abstandsbedingungen, simuliert mit XPBD.
//!
//! - **XPBD** (*extended position-based dynamics*, Macklin, Müller, Chentanez
//!   2016): Jedes Teilchen wird zuerst frei bewegt, dann werden die
//!   Bedingungen – hier feste Abstände zwischen Nachbarn – der Reihe nach
//!   durch Verschieben der Lagen erfüllt. Die Nachgiebigkeit
//!   (*compliance*) macht die Steifigkeit unabhängig von Schrittweite und
//!   Zahl der Durchläufe.
//! - **Hindernisse**: Ebene und Kugel als [`Collider`], beliebige ruhende
//!   Dreiecke über [`Cloth::add_obstacle_triangle`].
//! - **Abstoßung**: Kommt ein Teilchen einem fremden Dreieck näher als die
//!   Dicke des Tuchs, werden beide auseinandergeschoben.
//! - **Sicherung gegen Durchdringung**: Am Ende jedes Teilschritts prüft die
//!   kontinuierliche Kollisionserkennung ([`crate::ccd`]) jede Bewegung
//!   Punkt gegen Dreieck und Kante gegen Kante. Würden sich Teile
//!   durchdringen, wird die Bewegung der beteiligten Teilchen gekürzt: Sie
//!   halten kurz vor der Berührung an oder kehren auf ihre Lage vom Anfang
//!   des Teilschritts zurück; notfalls alle. War das Tuch vorher frei von
//!   Durchdringungen, ist es das danach auch – gleich wie schnell es sich
//!   bewegt und wie dünn es ist.
//!
//! Die Sicherung gilt für Dreiecke: das Tuch selbst und Hindernisse aus
//! Dreiecken. Ebene und Kugel schieben Teilchen nur hinaus.

mod collide;
mod diagnostics;
mod model;
mod step;

#[cfg(all(test, any(feature = "std", feature = "libm")))]
mod tests;

pub use model::{Cloth, ClothParams, Collider, StepStats, StepPhase};
pub use diagnostics::ClothDiagnostics;
