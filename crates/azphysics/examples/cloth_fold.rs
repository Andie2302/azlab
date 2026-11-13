//! Faltenwurf: Ein hochkant fallendes Tuch trifft auf den Boden und legt
//! sich in Falten auf sich selbst.
//!
//! Zeigt die Selbstkollision von [`azphysics::cloth::Cloth`]: Abstoßung
//! hält die Lagen auf Abstand, die kontinuierliche Kollisionserkennung
//! verhindert, dass sie sich durchdringen. Das Fenster zählt laufend die
//! Durchdringungen; die Zahl bleibt null.
//!
//! ```text
//! cargo demo cloth_fold
//! ```

#[path = "common/cloth_profile.rs"]
mod cloth_profile;

use azmath::Vector;
use azphysics::cloth::{Cloth, Collider};
use azviz::{Batch, Camera, Canvas, Color, Options, Scene, VizError};

const COLUMNS: usize = 28;
const ROWS: usize = 40;

struct Fold {
    camera: Camera,
    profiler: cloth_profile::ClothProfiler,
    show_ground: bool,
    cloth: Cloth<f32>,
    ground: [Collider<f32>; 1],
}

impl Fold {
    fn new() -> Self {
        let spacing = 0.035;
        // Fast senkrecht, leicht geneigt: So kippt das Tuch beim Auftreffen
        // und faltet sich.
        let across = Vector([spacing, 0.0, 0.0]);
        let down = Vector([0.0, -0.97 * spacing, 0.24 * spacing]);
        let origin = Vector([-0.5 * spacing * (COLUMNS - 1) as f32, 1.75, -0.15]);
        let mut cloth = Cloth::grid(COLUMNS, ROWS, origin, across, down, 0.004);
        cloth.params.thickness = 0.012;
        cloth.params.bend_compliance = 2.0e-2;
        cloth.params.friction = 0.6;
        let ground = [Collider::Plane { point: Vector([0.0, 0.0, 0.0]), normal: Vector([0.0, 1.0, 0.0]) }];
        Self { profiler: cloth_profile::ClothProfiler::default(), cloth, ground, camera: Camera::look_at([1.6, 1.2, 2.0], [0.0, 0.55, 0.15], 40.0), show_ground: true }
    }
}

impl Scene for Fold {
    fn title(&self) -> String {
        String::from("azphysics – Faltenwurf")
    }

    fn update(&mut self, dt: f32) {
        self.profiler.step(&mut self.cloth, dt, &self.ground);
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
        let positions: Vec<[f32; 3]> = self.cloth.positions().iter().map(|point| point.0).collect();
        batch.mesh(&positions, self.cloth.triangles(), Color::BLUE, Color::GREEN);
        batch.draw(canvas, &camera, true);
        let count = format!("Durchdringungen: {}", self.cloth.intersections());
        canvas.text([12.0, 46.0], &count, 20.0, Color::WHITE);
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(Fold::new(), Options { frames: 210, ..Options::default() })
}
