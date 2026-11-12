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

use azmath::Vector;
use azphysics::cloth::{Cloth, Collider};
use azviz::{Batch, Camera, Canvas, Color, Options, Scene, VizError};

const SIZE: usize = 34;
const RADIUS: f32 = 0.45;

struct Drape {
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
        Self { cloth, colliders }
    }
}

impl Scene for Drape {
    fn title(&self) -> String {
        String::from("azphysics – Tuch über Kugel")
    }

    fn update(&mut self, dt: f32) {
        self.cloth.step(dt, &self.colliders);
    }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        let camera = Camera::look_at([2.3, 1.9, 3.0], [0.0, 0.45, 0.0], 40.0);
        let mut batch = Batch::new();
        batch.floor(0.0, 0.5, 10, Color::rgb(60, 66, 80), Color::rgb(48, 53, 66));
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
