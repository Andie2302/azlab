//! Zeichenfläche im Speicher, für die Ausgabe in Dateien.

use crate::canvas::Canvas;
use crate::color::Color;
use tiny_skia::{FillRule, LineCap, Paint, PathBuilder, Pixmap, Rect, Stroke, Transform};

/// Zeichnet mit Kantenglättung in ein Bild im Speicher.
pub(crate) struct RasterCanvas {
    pixmap: Pixmap,
}

impl RasterCanvas {
    /// Neue Fläche; `None` bei Größe null.
    pub(crate) fn new(width: u32, height: u32) -> Option<Self> {
        Pixmap::new(width, height).map(|pixmap| Self { pixmap })
    }

    /// Das Bild.
    pub(crate) fn pixmap(&self) -> &Pixmap {
        &self.pixmap
    }

    fn paint(color: Color) -> Paint<'static> {
        let mut paint = Paint::default();
        paint.set_color_rgba8(color.r, color.g, color.b, color.a);
        paint.anti_alias = true;
        paint
    }
}

impl Canvas for RasterCanvas {
    fn width(&self) -> f32 {
        self.pixmap.width() as f32
    }

    fn height(&self) -> f32 {
        self.pixmap.height() as f32
    }

    fn clear(&mut self, color: Color) {
        self.pixmap.fill(tiny_skia::Color::from_rgba8(color.r, color.g, color.b, color.a));
    }

    fn line(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: Color) {
        let mut builder = PathBuilder::new();
        builder.move_to(from[0], from[1]);
        builder.line_to(to[0], to[1]);
        if let Some(path) = builder.finish() {
            let stroke = Stroke { width, line_cap: LineCap::Round, ..Stroke::default() };
            self.pixmap.stroke_path(&path, &Self::paint(color), &stroke, Transform::identity(), None);
        }
    }

    fn circle(&mut self, center: [f32; 2], radius: f32, color: Color) {
        if let Some(path) = PathBuilder::from_circle(center[0], center[1], radius) {
            self.pixmap.fill_path(&path, &Self::paint(color), FillRule::Winding, Transform::identity(), None);
        }
    }

    fn triangle(&mut self, a: [f32; 2], b: [f32; 2], c: [f32; 2], color: Color) {
        let mut builder = PathBuilder::new();
        builder.move_to(a[0], a[1]);
        builder.line_to(b[0], b[1]);
        builder.line_to(c[0], c[1]);
        builder.close();
        if let Some(path) = builder.finish() {
            self.pixmap.fill_path(&path, &Self::paint(color), FillRule::Winding, Transform::identity(), None);
        }
    }

    fn rect(&mut self, corner: [f32; 2], size: [f32; 2], color: Color) {
        if let Some(rect) = Rect::from_xywh(corner[0], corner[1], size[0], size[1]) {
            let mut paint = Self::paint(color);
            // Kacheln sollen lückenlos aneinanderstoßen.
            paint.anti_alias = false;
            self.pixmap.fill_rect(rect, &paint, Transform::identity(), None);
        }
    }
}
