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
//!
//! Die Geometrie und die Kollisionserkennung laufen in `no_std` ohne Heap.
//! Die Simulationen verwalten viele Teilchen und stehen deshalb hinter dem
//! Feature `alloc`.

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
pub mod geometry;
