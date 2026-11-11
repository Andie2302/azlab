//Regel #1 alles was nicht hinter feature flags steht MUSS in no_std ohne features laufen können.
//Der gesamte Workspace ist im NEUEN RUST STYLE ZU SCHREIBEN. Keine mod.rs irgendwo, im lib.rs darf KEIN CODE stehen!
//Bevor auch nur an Code Erstellung gedacht wird, muss im todo ordner alles geplant werden. inklusive abhängigkeiten, feature flags, zielen und nicht zielen. (DAS GILT GANZ BESONDERS FÜR KI!!!!)
//KI FÜHRT VON SICH AUS KEINEN CODE AUS! NICHT ZUM TESTEN, NICHT ZUM RECHNEN: NIEMALS!!! AUCH NICHT ZUM SCHNELLEN PRÜFEN!


//! # azmath
//!
//! Reine Mathematik-Bibliothek des azlab-Workspace: Zahlentypen,
//! Rechenmethoden und Löser. Fachliche Crates (Stoffdaten, Strömung,
//! Robotik, KI, ...) bauen darauf auf, azmath selbst kennt keine Fachdomäne.
//!
//! ## Grundsätze
//!
//! - `no_std` ohne `alloc` ist der Normalfall. Alles, was Heap oder `std`
//!   braucht, liegt hinter den Features `alloc` bzw. `std`.
//! - Algorithmen werden generisch über die Traits aus [`scalar`] geschrieben,
//!   nie direkt gegen `f32`/`f64`. Neue Zahlentypen (Dualzahlen, komplexe
//!   Zahlen, später `f16`/`f128`) funktionieren dadurch ohne Änderung an
//!   bestehendem Code.
//! - Wer die Elementarfunktionen tatsächlich berechnet (`libm`, `std`, ...),
//!   ist ein Implementierungsdetail des privaten Backend-Moduls und wird nur
//!   über Feature-Flags gewählt.
//!
//! ## Features
//!
//! | Feature | Wirkung |
//! |---------|---------|
//! | `libm` (Standard) | Elementarfunktionen über die Crate `libm`, läuft in `no_std` |
//! | `std`   | Elementarfunktionen über `std`, aktiviert auch `alloc` |
//! | `alloc` | Heap-basierte Typen (dynamische Größen) |

#![no_std]

#[cfg(feature = "std")]
extern crate std;

mod backend;
pub mod scalar;

pub use scalar::{Elementary, Field, Real, RealElementary, Ring};
