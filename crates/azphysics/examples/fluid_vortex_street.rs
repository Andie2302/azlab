//! Wirbelstraße: Ein Zylinder im Kanal wird von links angeströmt. Hinter ihm
//! lösen sich abwechselnd oben und unten Wirbel ab – die Kármánsche
//! Wirbelstraße.
//!
//! Zeigt [`azphysics::fluid::GridFluid`] mit Einlass, offenem Auslass,
//! Hindernis und Zähigkeit. Die Farbe ist die Wirbelstärke: Rot dreht gegen
//! den Uhrzeigersinn, Blau im Uhrzeigersinn.
//!
//! ```text
//! cargo demo fluid_vortex_street
//! ```

use azmath::Vector;
use azphysics::fluid::{Edges, GridFluid};
use azviz::{Canvas, Color, Options, Plane, Scene, VizError};

const NX: usize = 220;
const NY: usize = 88;
const CELL: f32 = 1.0 / NY as f32;
const SPEED: f32 = 1.0;
const RADIUS: f32 = 0.06;

struct VortexStreet {
    fluid: GridFluid<f32>,
}

impl VortexStreet {
    fn new() -> Self {
        let mut fluid = GridFluid::new(NX, NY, CELL, Edges::channel(SPEED));
        // Reynolds-Zahl = Geschwindigkeit · Durchmesser / Zähigkeit ≈ 240.
        fluid.params.viscosity = 5.0e-4;
        // Ein wenig aus der Mitte gerückt: Das stößt die Ablösung an.
        fluid.add_solid_disc(Vector([0.45, 0.515]), RADIUS);
        Self { fluid }
    }
}

impl Scene for VortexStreet {
    fn title(&self) -> String {
        String::from("azphysics – Kármánsche Wirbelstraße")
    }

    fn update(&mut self, dt: f32) {
        self.fluid.step(dt);
    }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        let plane = Plane::fit(canvas, [0.0, 0.0], [NX as f32 * CELL, NY as f32 * CELL], 16.0);
        let (calm, left, right) = (Color::rgb(24, 28, 38), Color::rgb(240, 90, 70), Color::rgb(70, 140, 250));
        let strongest = 6.0 * SPEED / RADIUS / 4.0;
        for j in 0..NY {
            for i in 0..NX {
                let color = if self.fluid.is_solid(i, j) {
                    Color::WHITE
                } else {
                    let spin = self.fluid.vorticity(i, j) / strongest;
                    calm.mix(if spin > 0.0 { left } else { right }, spin.abs().sqrt())
                };
                let min = [i as f32 * CELL, j as f32 * CELL];
                plane.rect(canvas, min, [min[0] + CELL, min[1] + CELL], color);
            }
        }
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(VortexStreet::new(), Options { width: 1000, height: 420, frames: 900, ..Options::default() })
}
