//! Zeichenfläche im Fenster.

use crate::canvas::Canvas;
use crate::color::Color;
use macroquad::prelude as mq;

/// Zeichnet unmittelbar in das offene Fenster.
pub(crate) struct WindowCanvas;

fn convert(color: Color) -> mq::Color {
    mq::Color::from_rgba(color.r, color.g, color.b, color.a)
}

impl Canvas for WindowCanvas {
    fn width(&self) -> f32 {
        mq::screen_width()
    }

    fn height(&self) -> f32 {
        mq::screen_height()
    }

    fn clear(&mut self, color: Color) {
        mq::clear_background(convert(color));
    }

    fn line(&mut self, from: [f32; 2], to: [f32; 2], width: f32, color: Color) {
        mq::draw_line(from[0], from[1], to[0], to[1], width, convert(color));
    }

    fn circle(&mut self, center: [f32; 2], radius: f32, color: Color) {
        mq::draw_circle(center[0], center[1], radius, convert(color));
    }

    fn triangle(&mut self, a: [f32; 2], b: [f32; 2], c: [f32; 2], color: Color) {
        mq::draw_triangle(
            mq::Vec2::new(a[0], a[1]),
            mq::Vec2::new(b[0], b[1]),
            mq::Vec2::new(c[0], c[1]),
            convert(color),
        );
    }

    fn rect(&mut self, corner: [f32; 2], size: [f32; 2], color: Color) {
        mq::draw_rectangle(corner[0], corner[1], size[0], size[1], convert(color));
    }

    fn text(&mut self, position: [f32; 2], text: &str, size: f32, color: Color) {
        mq::draw_text(text, position[0], position[1], size, convert(color));
    }
}
