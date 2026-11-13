//! Zeichenfläche im Fenster.

use crate::canvas::Canvas;
use crate::color::Color;
use macroquad::prelude as mq;

/// Zeichnet unmittelbar in das offene Fenster.
#[derive(Default)]
pub(crate) struct WindowCanvas {
    depth: Option<crate::DepthImage>,
    textures: alloc::vec::Vec<(u16, u16, mq::Texture2D)>,
    image_index: usize,
}

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
        self.image_index = 0;
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

    fn acquire_depth_image(&mut self, width: u16, height: u16) -> crate::DepthImage {
        let mut image = self.depth.take().unwrap_or_else(|| crate::DepthImage::new(0, 0));
        image.prepare(usize::from(width), usize::from(height));
        image
    }

    fn present_depth_image(&mut self, image: crate::DepthImage) {
        self.rgba_image(image.width as u16, image.height as u16, &image.rgba);
        self.depth = Some(image);
    }

    fn rgba_image(&mut self, width: u16, height: u16, rgba: &[u8]) {
        if rgba.len() != usize::from(width) * usize::from(height) * 4 { return; }
        // Jede Bildübertragung eines Frames erhält einen eigenen Slot:
        // ein späterer Upload darf vorher eingereihte Draws nicht verändern.
        let slot = self.image_index;
        self.image_index += 1;
        match self.textures.get(slot) {
            Some((w, h, texture)) if (*w, *h) == (width, height) => {
                texture.update_from_bytes(u32::from(width), u32::from(height), rgba);
            }
            _ => {
                let texture = mq::Texture2D::from_rgba8(width, height, rgba);
                texture.set_filter(mq::FilterMode::Nearest);
                if slot == self.textures.len() {
                    self.textures.push((width, height, texture));
                } else {
                    self.textures[slot] = (width, height, texture);
                }
            }
        }
        if let Some((_, _, texture)) = self.textures.get(slot) {
            mq::draw_texture(texture, 0.0, 0.0, mq::WHITE);
        }
    }

    fn text(&mut self, position: [f32; 2], text: &str, size: f32, color: Color) {
        mq::draw_text(text, position[0], position[1], size, convert(color));
    }
}
