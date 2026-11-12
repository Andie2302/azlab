//! Gemeinsamer Software-Tiefenpuffer für Fenster und Dateiausgabe.
use alloc::vec;
use alloc::vec::Vec;
use crate::Color;

pub(crate) struct DepthImage {
    pub(crate) width: usize,
    pub(crate) height: usize,
    pub(crate) rgba: Vec<u8>,
    depth: Vec<f32>,
}

impl DepthImage {
    pub(crate) fn new(width: usize, height: usize) -> Self {
        Self { width, height, rgba: vec![0; width * height * 4],
            depth: vec![f32::INFINITY; width * height] }
    }

    /// Bildkoordinaten und positive Kameratiefe je Ecke.
    pub(crate) fn triangle(&mut self, vertices: [[f32; 3]; 3], color: Color, priority: f32, wire: bool) {
        if color.a == 0 || !priority.is_finite()
            || vertices.iter().flatten().any(|v| !v.is_finite())
            || vertices.iter().any(|v| v[2] <= 0.0) {
            return;
        }
        let [a, b, c] = vertices;
        let edge = |p: [f32; 3], q: [f32; 3], x: f32, y: f32|
            (q[0] - p[0]) * (y - p[1]) - (q[1] - p[1]) * (x - p[0]);
        let area = edge(a, b, c[0], c[1]);
        if !area.is_finite() || area.abs() <= f32::MIN_POSITIVE { return; }
        let bounds = |axis: usize, limit: usize| {
            let lo = a[axis].min(b[axis]).min(c[axis]).floor().max(0.0).min(limit as f32) as usize;
            let hi = a[axis].max(b[axis]).max(c[axis]).ceil().max(0.0).min(limit as f32) as usize;
            (lo, hi)
        };
        let (left, right) = bounds(0, self.width);
        let (top, bottom) = bounds(1, self.height);
        let lengths = [[b, c], [c, a], [a, b]].map(|[p, q]| {
            ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2)).sqrt()
        });
        let dark = color.shade(0.75);
        for y in top..bottom {
            for x in left..right {
                let (px, py) = (x as f32 + 0.5, y as f32 + 0.5);
                let weights = [edge(b, c, px, py) / area,
                    edge(c, a, px, py) / area, edge(a, b, px, py) / area];
                if weights.iter().any(|&w| w < 0.0) { continue; }
                // 1/z ist nach der Perspektivprojektion linear im Bild.
                let inverse = weights[0] / a[2] + weights[1] / b[2] + weights[2] / c[2];
                let z = inverse.recip() - priority;
                let index = y * self.width + x;
                let outline = wire && (0..3).any(|i|
                    weights[i] * area.abs() <= 0.5 * lengths[i]);
                let pixel = if outline { dark } else { color };
                let rgba = [pixel.r, pixel.g, pixel.b, pixel.a];
                let old = &self.rgba[index * 4..index * 4 + 4];
                // Exakt gleiche Tiefe: deterministische Wahl statt Zeichenreihenfolge.
                if !z.is_finite() || z > self.depth[index]
                    || (z == self.depth[index] && old <= rgba.as_slice()) { continue; }
                self.depth[index] = z;
                self.rgba[index * 4..index * 4 + 4]
                    .copy_from_slice(&[pixel.r, pixel.g, pixel.b, pixel.a]);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crossing_surfaces_are_order_independent() {
        // Gleiche Projektion, aber wechselnde Tiefenreihenfolge.
        let a = [[0.0, 0.0, 1.0], [8.0, 0.0, 9.0], [0.0, 8.0, 1.0]];
        let b = [[0.0, 0.0, 2.0], [8.0, 0.0, 2.0], [0.0, 8.0, 2.0]];
        let render = |reverse: bool, wire: bool| {
            let mut image = DepthImage::new(8, 8);
            let surfaces = if reverse { [(b, Color::BLUE), (a, Color::RED)] }
                else { [(a, Color::RED), (b, Color::BLUE)] };
            for (tri, color) in surfaces { image.triangle(tri, color, 0.0, wire); }
            image
        };
        for wire in [false, true] {
            assert_eq!(render(false, wire).rgba, render(true, wire).rgba);
        }
        let image = render(false, false);
        assert_eq!(&image.rgba[0..4], &[Color::RED.r, Color::RED.g, Color::RED.b, 255]);
        assert_eq!(&image.rgba[6 * 4..7 * 4], &[Color::BLUE.r, Color::BLUE.g, Color::BLUE.b, 255]);
        // Perspektivisch korrekt: am linken Pixel ist z etwa 1.06, nicht 1.5.
        assert!((image.depth[0] - 1.0588236).abs() < 1.0e-5);
    }

    #[test]
    fn clipping_degenerate_geometry_and_clear() {
        let mut image = DepthImage::new(4, 4);
        image.triangle([[-20.0, -20.0, 2.0], [20.0, -20.0, 2.0], [0.0, 20.0, 2.0]],
            Color::GREEN, 0.0, false);
        assert!(image.rgba.chunks_exact(4).all(|p| p[3] == 255));
        let before = image.rgba.clone();
        image.triangle([[1.0, 1.0, 1.0]; 3], Color::RED, 0.0, true);
        image.triangle([[f32::NAN, 1.0, 1.0]; 3], Color::RED, 0.0, false);
        assert_eq!(image.rgba, before);
        assert!(DepthImage::new(4, 4).rgba.iter().all(|&v| v == 0));
    }
}
