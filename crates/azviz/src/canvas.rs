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

    /// RGBA-Bild ohne vormultipliziertes Alpha, ab der linken oberen Ecke.
    /// Standardweg für eigene Zeichenflächen: gleichfarbige Zeilenabschnitte.
    fn rgba_image(&mut self, width: u16, height: u16, rgba: &[u8]) {
        let width = usize::from(width);
        if rgba.len() != width * usize::from(height) * 4 { return; }
        for y in 0..usize::from(height) {
            let mut x = 0;
            while x < width {
                let at = (y * width + x) * 4;
                let p = &rgba[at..at + 4];
                let mut end = x + 1;
                while end < width && rgba[(y * width + end) * 4..(y * width + end + 1) * 4] == *p {
                    end += 1;
                }
                if p[3] != 0 {
                    self.rect([x as f32, y as f32], [(end - x) as f32, 1.0],
                        Color { r: p[0], g: p[1], b: p[2], a: p[3] });
                }
                x = end;
            }
        }
    }

    /// Text; nur das Fenster stellt ihn dar, Dateien lassen ihn weg.
    fn text(&mut self, _position: [f32; 2], _text: &str, _size: f32, _color: Color) {}
}
