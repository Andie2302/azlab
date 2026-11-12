//! Rauch: Eine warme Rauchfahne steigt in einem geschlossenen Kasten auf und
//! umströmt ein Hindernis.
//!
//! Zeigt [`azphysics::fluid::GridFluid`]: Auftrieb, Transport und die
//! Druckkorrektur, die die Strömung quellenfrei hält.
//!
//! ```text
//! cargo demo fluid_smoke
//! ```

use azmath::Vector;
use azphysics::fluid::{Edges, GridFluid};
use azviz::{Canvas, Color, Options, Plane, Scene, VizError};

const NX: usize = 120;
const NY: usize = 160;
const CELL: f32 = 1.0 / NX as f32;

struct Smoke {
    fluid: GridFluid<f32>,
    time: f32,
}

impl Smoke {
    fn new() -> Self {
        let mut fluid = GridFluid::new(NX, NY, CELL, Edges::closed());
        fluid.params.buoyancy = 3.0;
        fluid.params.dissipation = 0.15;
        fluid.add_solid_disc(Vector([0.56, 0.62]), 0.09);
        Self { fluid, time: 0.0 }
    }
}

impl Scene for Smoke {
    fn title(&self) -> String {
        String::from("azphysics – Rauch (Navier-Stokes auf dem Gitter)")
    }

    fn update(&mut self, dt: f32) {
        self.time += dt;
        // Quelle am Boden; sie pendelt leicht, damit die Fahne lebendig wird.
        let centre = NX as f32 * (0.5 + 0.03 * (self.time * 1.7).sin());
        for j in 2..8 {
            for i in 0..NX {
                let distance = (i as f32 + 0.5 - centre).abs() / (NX as f32 * 0.05);
                if distance < 1.0 {
                    self.fluid.add_smoke(i, j, (1.0 - self.fluid.smoke_at(i, j)).max(0.0) * 0.5);
                    self.fluid.add_velocity(i, j, Vector([0.0, 0.6 * dt]));
                }
            }
        }
        self.fluid.step(dt);
    }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        let plane = Plane::fit(canvas, [0.0, 0.0], [NX as f32 * CELL, NY as f32 * CELL], 16.0);
        let (thin, dense) = (Color::rgb(24, 28, 38), Color::rgb(250, 236, 214));
        for j in 0..NY {
            for i in 0..NX {
                let color = if self.fluid.is_solid(i, j) {
                    Color::GRAY
                } else {
                    thin.mix(dense, self.fluid.smoke_at(i, j).sqrt())
                };
                let min = [i as f32 * CELL, j as f32 * CELL];
                plane.rect(canvas, min, [min[0] + CELL, min[1] + CELL], color);
            }
        }
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(Smoke::new(), Options { width: 600, height: 800, frames: 360, ..Options::default() })
}
