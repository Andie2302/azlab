//! Dammbruch: Eine Wassersäule aus Teilchen stürzt in ein Becken und
//! schwappt über ein Hindernis.
//!
//! Zeigt [`azphysics::fluid::ParticleFluid`]. Die Farbe gibt die
//! Geschwindigkeit an: dunkelblau ruht, weiß ist schnell.
//!
//! ```text
//! cargo demo fluid_dam_break
//! ```

use azmath::Vector;
use azphysics::fluid::ParticleFluid;
use azphysics::geometry::Aabb;
use azviz::{Canvas, Color, Options, Plane, Scene, VizError};

const WIDTH: f32 = 4.0;
const HEIGHT: f32 = 2.4;

struct DamBreak {
    fluid: ParticleFluid<f32>,
}

impl DamBreak {
    fn new() -> Self {
        let mut fluid = ParticleFluid::new(Aabb { min: Vector([0.0, 0.0]), max: Vector([WIDTH, HEIGHT]) }, 0.04);
        fluid.add_obstacle(Vector([2.2, 0.0]), Vector([2.4, 0.4]));
        fluid.add_block(Vector([0.0, 0.0]), Vector([1.2, 2.0]));
        Self { fluid }
    }
}

impl Scene for DamBreak {
    fn title(&self) -> String {
        format!("azphysics – Dammbruch ({} Teilchen)", self.fluid.positions().len())
    }

    fn update(&mut self, dt: f32) {
        self.fluid.step(dt);
    }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        let plane = Plane::fit(canvas, [0.0, 0.0], [WIDTH, HEIGHT], 24.0);
        // Becken und Hindernis.
        let wall = Color::GRAY;
        let corners = [[0.0, HEIGHT], [0.0, 0.0], [WIDTH, 0.0], [WIDTH, HEIGHT]];
        for pair in corners.windows(2) {
            canvas.line(plane.point(pair[0]), plane.point(pair[1]), 3.0, wall);
        }
        for block in self.fluid.obstacles() {
            plane.rect(canvas, block.min.0, block.max.0, wall);
        }
        let radius = plane.length(self.fluid.spacing() * 0.55);
        let (slow, fast) = (Color::rgb(30, 90, 200), Color::rgb(235, 245, 255));
        for (point, velocity) in self.fluid.positions().iter().zip(self.fluid.velocities()) {
            let color = slow.mix(fast, velocity.norm() / 5.0);
            canvas.circle(plane.point(point.0), radius, color);
        }
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(DamBreak::new(), Options { width: 900, height: 560, ..Options::default() })
}
