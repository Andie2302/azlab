//! Tuch über Kugel: Ein Tuch fällt auf eine Kugel, legt sich in Falten
//! darüber und rutscht langsam ab.
//!
//! Zeigt [`azphysics::cloth::Cloth`]: XPBD-Abstandsbedingungen, Hindernisse
//! mit Reibung und die Selbstkollision, die verhindert, dass sich die Falten
//! durchdringen.
//!
//! ```text
//! cargo demo cloth_drape
//! ```

#[path = "common/cloth_profile.rs"]
mod cloth_profile;

use azmath::Vector;
use azphysics::cloth::{Cloth, Collider};
use azviz::{Batch, Camera, Canvas, Color, Options, Scene, VizError};

const SIZE: usize = 34;
const RADIUS: f32 = 0.45;

struct Drape {
    camera: Camera,
    profiler: cloth_profile::ClothProfiler,
    show_ground: bool,
    cloth: Cloth<f32>,
    colliders: [Collider<f32>; 2],
}

impl Drape {
    fn new() -> Self {
        let spacing = 1.8 / (SIZE - 1) as f32;
        // Leicht verdreht und aus der Mitte gerückt, damit der Fall nicht
        // spiegelbildlich bleibt.
        let (across, down) = (Vector([spacing, 0.0, 0.25 * spacing]), Vector([-0.25 * spacing, 0.0, spacing]));
        let mut cloth = Cloth::grid(SIZE, SIZE, Vector([-0.62, 1.35, -1.0]), across, down, 0.005);
        cloth.params.thickness = 0.012;
        //cloth.params.gravity = Vector([0.0, 0.0, 0.0]);
        cloth.params.bend_compliance = 5.0e-3;
        let colliders = [
            Collider::Sphere { center: Vector([0.0, RADIUS, 0.0]), radius: RADIUS },
            Collider::Plane { point: Vector([0.0, 0.0, 0.0]), normal: Vector([0.0, 1.0, 0.0]) },
        ];
        Self { profiler: cloth_profile::ClothProfiler::default(), cloth, colliders, camera: Camera::look_at([2.3, 1.9, 3.0], [0.0, 0.45, 0.0], 40.0), show_ground: true }
    }
}

impl Scene for Drape {
    fn title(&self) -> String {
        String::from("azphysics – Tuch über Kugel")
    }

    fn update(&mut self, dt: f32) {
        self.profiler.step(&mut self.cloth, dt, &self.colliders);
    }

    fn set_profiling(&mut self, enabled: bool) { self.profiler.enabled = enabled; }

    fn profile_details(&mut self) -> Option<String> { self.profiler.report() }

    fn orbit(&mut self, yaw: f32, pitch: f32, zoom: f32) {
        self.camera.orbit(yaw, pitch, zoom);
    }

    fn toggle_ground(&mut self) { self.show_ground = !self.show_ground; }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        let camera = self.camera;
        let mut batch = Batch::new();
        if self.show_ground { batch.floor(0.0, 0.5, 10, Color::rgb(60, 66, 80), Color::rgb(48, 53, 66)); }
        // Etwas kleiner als das Hindernis, damit das Tuch sichtbar aufliegt.
        batch.sphere([0.0, RADIUS, 0.0], RADIUS * 0.97, Color::GRAY);
        let positions: Vec<[f32; 3]> = self.cloth.positions().iter().map(|point| point.0).collect();
        batch.mesh(&positions, self.cloth.triangles(), Color::ORANGE, Color::RED);
        batch.draw(canvas, &camera, true);
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(Drape::new(), Options { frames: 240, ..Options::default() })
}
