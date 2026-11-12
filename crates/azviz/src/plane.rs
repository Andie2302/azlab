//! Blick auf die Ebene: Abbildung eines Rechtecks der Welt auf die
//! Zeichenfläche.

use crate::canvas::Canvas;

/// Bildet einen rechteckigen Ausschnitt der Ebene (`y` nach oben) so auf die
/// Zeichenfläche ab, dass er ganz hineinpasst und nicht verzerrt wird.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plane {
    /// Bildpunkt der linken unteren Ecke des Ausschnitts.
    corner: [f32; 2],
    /// Bildpunkte je Längeneinheit.
    scale: f32,
    min: [f32; 2],
}

impl Plane {
    /// Passt den Ausschnitt von `min` bis `max` mit `margin` Bildpunkten
    /// Rand in die Zeichenfläche ein.
    pub fn fit(canvas: &dyn Canvas, min: [f32; 2], max: [f32; 2], margin: f32) -> Self {
        let (width, height) = ((max[0] - min[0]).max(f32::MIN_POSITIVE), (max[1] - min[1]).max(f32::MIN_POSITIVE));
        let scale = ((canvas.width() - 2.0 * margin) / width).min((canvas.height() - 2.0 * margin) / height);
        let corner = [(canvas.width() - width * scale) * 0.5, (canvas.height() + height * scale) * 0.5];
        Self { corner, scale, min }
    }

    /// Bildpunkt eines Punkts der Ebene.
    pub fn point(&self, point: [f32; 2]) -> [f32; 2] {
        [self.corner[0] + (point[0] - self.min[0]) * self.scale, self.corner[1] - (point[1] - self.min[1]) * self.scale]
    }

    /// Länge in Bildpunkten.
    pub fn length(&self, length: f32) -> f32 {
        length * self.scale
    }

    /// Zeichnet ein achsenparalleles Rechteck der Ebene.
    pub fn rect(&self, canvas: &mut dyn Canvas, min: [f32; 2], max: [f32; 2], color: crate::Color) {
        let (top_left, bottom_right) = (self.point([min[0], max[1]]), self.point([max[0], min[1]]));
        // Auf ganze Bildpunkte nach außen runden, damit zwischen Kacheln
        // keine Fugen bleiben.
        let top_left = [top_left[0].floor(), top_left[1].floor()];
        let size = [bottom_right[0].ceil() - top_left[0], bottom_right[1].ceil() - top_left[1]];
        canvas.rect(top_left, size, color);
    }
}
