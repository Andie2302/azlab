//! Die Zeichenfläche.

use crate::color::Color;

/// Zeichenfläche in Bildpunkten; der Ursprung liegt oben links, `y` wächst
/// nach unten.
pub trait Canvas {
    /// Breite in Bildpunkten.
    fn width(&self) -> f32;

    /// Höhe in Bildpunkten.
    fn height(&self) -> f32;

    /// Füllt die ganze Fläche.
    fn clear(&mut self, color: Color);

    /// Strecke der Breite `width`.
    fn line(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: Color);

    /// Gefüllter Kreis.
    fn circle(&mut self, center: [f32; 2], radius: f32, color: Color);

    /// Gefülltes Dreieck.
    fn triangle(&mut self, a: [f32; 2], b: [f32; 2], c: [f32; 2], color: Color);

    /// Gefülltes achsenparalleles Rechteck ab der Ecke `corner`.
    fn rect(&mut self, corner: [f32; 2], size: [f32; 2], color: Color);

    /// Text; nur das Fenster stellt ihn dar, Dateien lassen ihn weg.
    fn text(&mut self, _position: [f32; 2], _text: &str, _size: f32, _color: Color) {}
}
