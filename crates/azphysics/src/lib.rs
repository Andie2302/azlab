//! # azphysics
//!
//! Physik-Simulationen auf der Grundlage von azmath.
//!
//! | Modul | Inhalt | Braucht |
//! |-------|--------|---------|
//! | [`geometry`] | nächste Punkte und Abstände zwischen Punkt, Strecke und Dreieck; Schnitt-Tests; Hüllquader | – |
//! | [`ccd`] | kontinuierliche Kollisionserkennung: Wann berühren sich bewegte Punkte, Kanten und Dreiecke? | – |
//! | [`broadphase`] | Suchstruktur für Paare, die sich überhaupt nahe kommen können | `alloc` |
//! | [`cloth`] | Gewebe mit XPBD, Kollision mit Hindernissen und mit sich selbst, ohne Durchdringung | `alloc` |
//! | [`fluid`] | Navier-Stokes auf dem Gitter (Rauch, Umströmung) und Flüssigkeit aus Teilchen | `alloc` |
//!
//! Die Geometrie und die Kollisionserkennung laufen in `no_std` ohne Heap.
//! Die Simulationen verwalten viele Teilchen und stehen deshalb hinter dem
//! Feature `alloc`.
//!
//! Die Beispiele unter `examples/` zeigen die Simulationen in einem
//! Fenster, etwa `cargo demo cloth_drape`.
//!
//! # Beispiel
//!
//! Ein Tuch fällt auf den Boden; es durchdringt dabei weder ihn noch sich
//! selbst.
//!
//! ```
//! # #[cfg(all(feature = "alloc", any(feature = "std", feature = "libm")))]
//! # {
//! use azmath::Vector;
//! use azphysics::cloth::{Cloth, Collider};
//!
//! // Tuch aus 12 × 12 Teilchen, 1 m × 1 m, in 30 cm Höhe.
//! let spacing = 1.0 / 11.0;
//! let mut cloth: Cloth<f64> = Cloth::grid(
//!     12,
//!     12,
//!     Vector([-0.5, 0.3, -0.5]),
//!     Vector([spacing, 0.0, 0.0]),
//!     Vector([0.0, 0.0, spacing]),
//!     0.01,
//! );
//! let ground = [Collider::Plane { point: Vector([0.0; 3]), normal: Vector([0.0, 1.0, 0.0]) }];
//!
//! for _ in 0..40 {
//!     cloth.step(1.0 / 60.0, &ground);
//! }
//! // Das Tuch liegt auf dem Boden: knapp darüber, nirgends darunter.
//! assert!(cloth.positions().iter().all(|point| point[1] > 0.0 && point[1] < 0.05));
//! assert_eq!(cloth.intersections(), 0);
//! # }
//! ```

#![no_std]

#[cfg(feature = "alloc")]
extern crate alloc;

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "alloc")]
pub mod broadphase;
pub mod ccd;
#[cfg(feature = "alloc")]
pub mod cloth;
#[cfg(feature = "alloc")]
pub mod fluid;
pub mod geometry;
