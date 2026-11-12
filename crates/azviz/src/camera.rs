//! Blick in den Raum: Projektion und schattierte Dreiecksnetze.

use crate::canvas::Canvas;
use crate::color::Color;
use alloc::vec::Vec;

type Point = [f32; 3];

fn sub(a: Point, b: Point) -> Point {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: Point, b: Point) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Point, b: Point) -> Point {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn normalize(a: Point) -> Point {
    let length = dot(a, a).sqrt();
    if length > 0.0 { [a[0] / length, a[1] / length, a[2] / length] } else { [0.0, 1.0, 0.0] }
}

/// Kamera mit Zentralprojektion; `y` zeigt im Raum nach oben.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera {
    eye: Point,
    right: Point,
    up: Point,
    forward: Point,
    /// Brennweite relativ zur Bildhöhe; größer heißt stärker vergrößert.
    focal: f32,
}

impl Camera {
    /// Kamera am Ort `eye` mit Blick auf `target`; `field_of_view` ist der
    /// senkrechte Öffnungswinkel in Grad.
    pub fn look_at(eye: Point, target: Point, field_of_view: f32) -> Self {
        let forward = normalize(sub(target, eye));
        let right = normalize(cross(forward, [0.0, 1.0, 0.0]));
        let up = cross(right, forward);
        let focal = 0.5 / (field_of_view.to_radians() * 0.5).tan();
        Self { eye, right, up, forward, focal }
    }

    /// Bildpunkt und Tiefe eines Raumpunkts; `None`, wenn er hinter der
    /// Kamera liegt.
    pub fn project(&self, point: Point, canvas: &dyn Canvas) -> Option<([f32; 2], f32)> {
        let relative = sub(point, self.eye);
        let depth = dot(relative, self.forward);
        if depth <= 1.0e-3 {
            return None;
        }
        let scale = self.focal * canvas.height() / depth;
        Some((
            [
                canvas.width() * 0.5 + dot(relative, self.right) * scale,
                canvas.height() * 0.5 - dot(relative, self.up) * scale,
            ],
            depth,
        ))
    }
}

/// Sammlung von Dreiecken im Raum, die gemeinsam gezeichnet werden:
/// schattiert nach ihrer Neigung zum Licht und von hinten nach vorn sortiert,
/// damit Vorderes Hinteres verdeckt.
#[derive(Clone, Debug, Default)]
pub struct Batch {
    triangles: Vec<([Point; 3], Color, Color, f32)>,
    /// Vorrang der Dreiecke, die als Nächstes hinzukommen.
    priority: f32,
    /// Boden: liegt unter allem und wird deshalb zuerst gezeichnet.
    ground: Vec<([Point; 4], Color)>,
}

impl Batch {
    /// Leere Sammlung.
    pub fn new() -> Self {
        Self::default()
    }

    /// Gibt allen danach hinzugefügten Dreiecken Vorrang: Sie werden so
    /// einsortiert, als lägen sie um `distance` näher an der Kamera. Damit
    /// bleibt etwas sichtbar, das dicht auf einer anderen Fläche aufliegt
    /// (ein Tuch auf einer Platte). Null stellt den Normalfall wieder her.
    pub fn set_priority(&mut self, distance: f32) {
        self.priority = distance;
    }

    /// Ein Dreieck. `front` färbt die Seite, auf der die Ecken gegen den
    /// Uhrzeigersinn laufen, `back` die Rückseite.
    pub fn triangle(&mut self, corners: [Point; 3], front: Color, back: Color) {
        self.triangles.push((corners, front, back, self.priority));
    }

    /// Ein Dreiecksnetz aus Punkten und Dreiecken über deren Nummern.
    pub fn mesh(&mut self, positions: &[Point], triangles: &[[usize; 3]], front: Color, back: Color) {
        for triangle in triangles {
            let corners = [positions.get(triangle[0]), positions.get(triangle[1]), positions.get(triangle[2])];
            if let [Some(&a), Some(&b), Some(&c)] = corners {
                self.triangle([a, b, c], front, back);
            }
        }
    }

    /// Ein ebenes Viereck mit den Ecken reihum, zerlegt in
    /// `divisions × divisions` Felder. Kleine Felder ordnen sich beim
    /// Zeichnen richtig zwischen andere Dreiecke ein; ein einziges großes
    /// Viereck verdeckt sonst leicht, was davor liegt.
    pub fn quad(&mut self, corners: [Point; 4], divisions: usize, front: Color, back: Color) {
        let divisions = divisions.max(1);
        let point = |column: usize, row: usize| {
            let (s, t) = (column as f32 / divisions as f32, row as f32 / divisions as f32);
            core::array::from_fn(|axis| {
                let near = corners[0][axis] + (corners[1][axis] - corners[0][axis]) * s;
                let far = corners[3][axis] + (corners[2][axis] - corners[3][axis]) * s;
                near + (far - near) * t
            })
        };
        for row in 0..divisions {
            for column in 0..divisions {
                let (a, b) = (point(column, row), point(column + 1, row));
                let (c, d) = (point(column + 1, row + 1), point(column, row + 1));
                self.triangle([a, b, c], front, back);
                self.triangle([a, c, d], front, back);
            }
        }
    }

    /// Eine Kugel.
    pub fn sphere(&mut self, centre: Point, radius: f32, color: Color) {
        const RINGS: usize = 12;
        const SEGMENTS: usize = 24;
        let point = |ring: usize, segment: usize| {
            let polar = core::f32::consts::PI * ring as f32 / RINGS as f32;
            let around = core::f32::consts::TAU * segment as f32 / SEGMENTS as f32;
            [
                centre[0] + radius * polar.sin() * around.cos(),
                centre[1] + radius * polar.cos(),
                centre[2] + radius * polar.sin() * around.sin(),
            ]
        };
        for ring in 0..RINGS {
            for segment in 0..SEGMENTS {
                let (a, b) = (point(ring, segment), point(ring, segment + 1));
                let (c, d) = (point(ring + 1, segment), point(ring + 1, segment + 1));
                self.triangle([a, b, c], color, color);
                self.triangle([b, d, c], color, color);
            }
        }
    }

    /// Ein waagrechter Boden in der Höhe `height`, als Schachbrett aus
    /// `tiles × tiles` Feldern der Kantenlänge `tile` um den Ursprung. Er
    /// wird vor allem anderen gezeichnet, liegt also immer dahinter.
    pub fn floor(&mut self, height: f32, tile: f32, tiles: usize, light: Color, dark: Color) {
        let half = tiles as f32 * tile * 0.5;
        for row in 0..tiles {
            for column in 0..tiles {
                let (x, z) = (column as f32 * tile - half, row as f32 * tile - half);
                let color = if (row + column) % 2 == 0 { light } else { dark };
                let corners =
                    [[x, height, z], [x, height, z + tile], [x + tile, height, z + tile], [x + tile, height, z]];
                self.ground.push((corners, color));
            }
        }
    }

    /// Zeichnet alles; mit `wire` werden die Kanten nachgezogen.
    pub fn draw(&self, canvas: &mut dyn Canvas, camera: &Camera, wire: bool) {
        for &(corners, color) in &self.ground {
            let projected = corners.map(|corner| camera.project(corner, canvas));
            if let [Some(a), Some(b), Some(c), Some(d)] = projected {
                canvas.triangle(a.0, b.0, c.0, color);
                canvas.triangle(a.0, c.0, d.0, color);
                // Die Naht zwischen den beiden Dreiecken schließen.
                canvas.line(a.0, c.0, 1.5, color);
            }
        }
        let light = normalize([0.4, 1.0, 0.6]);
        let mut visible: Vec<(f32, [[f32; 2]; 3], Color)> = Vec::new();
        for &([a, b, c], front, back, priority) in &self.triangles {
            let (Some(pa), Some(pb), Some(pc)) =
                (camera.project(a, canvas), camera.project(b, canvas), camera.project(c, canvas))
            else {
                continue;
            };
            let normal = normalize(cross(sub(b, a), sub(c, a)));
            let facing = dot(normal, sub(camera.eye, a)) >= 0.0;
            // Beleuchtung der sichtbaren Seite.
            let lit = dot(if facing { normal } else { [-normal[0], -normal[1], -normal[2]] }, light);
            let brightness = 0.45 + 0.55 * lit.max(0.0);
            let color = if facing { front } else { back }.shade(brightness);
            visible.push(((pa.1 + pb.1 + pc.1) / 3.0 - priority, [pa.0, pb.0, pc.0], color));
        }
        visible.sort_by(|a, b| b.0.total_cmp(&a.0));
        for (_, [a, b, c], color) in visible {
            canvas.triangle(a, b, c, color);
            if wire {
                let edge = color.shade(0.75);
                canvas.line(a, b, 1.0, edge);
                canvas.line(b, c, 1.0, edge);
                canvas.line(c, a, 1.0, edge);
            }
        }
    }
}
