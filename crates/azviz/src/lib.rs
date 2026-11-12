//! # azviz
//!
//! Anzeige für die Beispiele des azlab-Workspace.
//!
//! Eine Simulation beschreibt sich als [`Scene`]: Sie rechnet sich
//! schrittweise weiter und zeichnet sich auf eine [`Canvas`]. [`run`] zeigt
//! sie dann entweder live in einem Fenster oder schreibt sie als Bild oder
//! Animation in eine Datei – je nach Kommandozeile:
//!
//! | Aufruf | Wirkung |
//! |--------|---------|
//! | *(ohne Argumente)* | Fenster; Leertaste hält an, Esc beendet |
//! | `--gif datei.gif` | Animation als GIF |
//! | `--png datei.png` | Standbild nach den gerechneten Bildern |
//! | `--frames N` | Zahl der Bilder für `--gif` und `--png` |
//! | `--size BxH` | Bildgröße in Bildpunkten, etwa `--size 480x360` |
//!
//! Die Ausgabe in Dateien braucht keinen Bildschirm.

extern crate alloc;

mod camera;
mod canvas;
mod color;
mod plane;
mod raster;
mod runner;
mod window;

#[cfg(test)]
mod tests;

pub use camera::{Batch, Camera};
pub use canvas::Canvas;
pub use color::Color;
pub use plane::Plane;
pub use runner::{Options, Scene, VizError, run};
