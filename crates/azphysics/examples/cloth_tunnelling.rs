//! Tunneln: Zwei gleiche Tücher werden mit 60 m/s auf eine dünne Platte
//! geschossen. Hinten (rot) ist die Sicherung gegen Durchdringung
//! ausgeschaltet – das Tuch überspringt die Platte zwischen zwei
//! Zeitschritten und fliegt hindurch. Vorn (grün) ist sie eingeschaltet –
//! das Tuch bleibt an der Platte hängen.
//!
//! Zeigt die kontinuierliche Kollisionserkennung ([`azphysics::ccd`]) im
//! Gewebe ([`azphysics::cloth::Cloth`]). Beide Tücher rechnen mit denselben
//! Einstellungen; nur `params.ccd` unterscheidet sie. Die Wiedergabe läuft
//! in zehnfacher Zeitlupe und beginnt alle vier Sekunden von vorn.
//!
//! ```text
//! cargo demo cloth_tunnelling
//! ```

#[path = "common/cloth_profile.rs"]
mod cloth_profile;

use azmath::Vector;
use azphysics::cloth::{Cloth, Collider};
use azviz::{Batch, Camera, Canvas, Color, Options, Scene, VizError};

const SIZE: usize = 20;
const SPEED: f32 = 60.0;
const SLOW_MOTION: f32 = 0.1;
/// Bilder bis zum Neustart.
const LOOP: u32 = 240;
/// Abstand der beiden Bahnen von der Mitte.
const LANE: f32 = 1.1;

struct Tunnelling {
    camera: Camera,
    profiler: cloth_profile::ClothProfiler,
    show_ground: bool,
    lanes: [Cloth<f32>; 2],
    ground: [Collider<f32>; 1],
    frame: u32,
}

/// Ecken der Platte einer Bahn bei `z`.
fn plate(z: f32) -> [[f32; 3]; 4] {
    [[0.0, 0.2, z - 0.7], [0.0, 0.2, z + 0.7], [0.0, 1.6, z + 0.7], [0.0, 1.6, z - 0.7]]
}

/// Tuch und Platte einer Bahn bei `z`.
fn lane(z: f32, ccd: bool) -> Cloth<f32> {
    let spacing = 0.8 / (SIZE - 1) as f32;
    // Senkrechtes Tuch, quer zur Flugrichtung.
    let mut cloth = Cloth::grid(
        SIZE,
        SIZE,
        Vector([-2.2, 1.3, z - 0.4]),
        Vector([0.0, 0.0, spacing]),
        Vector([0.0, -spacing, 0.0]),
        0.005,
    );
    cloth.params.ccd = ccd;
    // Ein Teilschritt je Zeitschritt: Das Tuch legt dazwischen zehn
    // Zentimeter zurück, das Zehnfache seiner Dicke.
    cloth.params.substeps = 1;
    cloth.params.iterations = 8;
    cloth.params.thickness = 0.01;
    cloth.set_all_velocities(Vector([SPEED, 0.0, 0.0]));
    // Die Platte: zwei Dreiecke in der Ebene x = 0.
    let corners = plate(z).map(Vector);
    cloth.add_obstacle_triangle(corners[0], corners[1], corners[2]);
    cloth.add_obstacle_triangle(corners[0], corners[2], corners[3]);
    cloth
}

impl Tunnelling {
    fn new() -> Self {
        Self {
            lanes: [lane(-LANE, false), lane(LANE, true)],
            ground: [Collider::Plane { point: Vector([0.0, 0.0, 0.0]), normal: Vector([0.0, 1.0, 0.0]) }],
            frame: 0,
            profiler: cloth_profile::ClothProfiler::default(),
            camera: Camera::look_at([-3.0, 2.3, 4.9], [-0.3, 0.75, 0.0], 36.0),
            show_ground: true,
        }
    }
}

impl Scene for Tunnelling {
    fn title(&self) -> String {
        String::from("azphysics – Tunneln: hinten ohne, vorn mit Sicherung")
    }

    fn update(&mut self, dt: f32) {
        self.frame += 1;
        if self.frame >= LOOP {
            let (camera, show_ground) = (self.camera, self.show_ground);
            let profiler = core::mem::take(&mut self.profiler);
            *self = Self::new();
            self.profiler = profiler;
            self.camera = camera;
            self.show_ground = show_ground;
        }
        for cloth in &mut self.lanes {
            self.profiler.step(cloth, dt * SLOW_MOTION, &self.ground);
        }
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
        if self.show_ground { batch.floor(0.0, 0.5, 16, Color::rgb(60, 66, 80), Color::rgb(48, 53, 66)); }
        let colors = [(Color::RED, Color::ORANGE), (Color::GREEN, Color::rgb(170, 220, 120))];
        for ((cloth, (front, back)), z) in self.lanes.iter().zip(colors).zip([-LANE, LANE]) {
            let positions: Vec<[f32; 3]> = cloth.positions().iter().map(|point| point.0).collect();
            batch.quad(plate(z), 14, Color::GRAY, Color::GRAY);
            // Das Tuch liegt dicht an der Platte an und soll davor sichtbar
            // bleiben.
            batch.set_priority(0.08);
            batch.mesh(&positions, &cloth.triangles()[..cloth.cloth_triangle_count()], front, back);
            batch.set_priority(0.0);
        }
        batch.draw(canvas, &camera, true);
        for (z, label) in [(-LANE, "ohne Sicherung"), (LANE, "mit Sicherung (CCD)")] {
            if let Some((point, _)) = camera.project([-2.2, 1.6, z], canvas) {
                canvas.text(point, label, 22.0, Color::WHITE);
            }
        }
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(Tunnelling::new(), Options { width: 960, height: 600, frames: LOOP, ..Options::default() })
}
