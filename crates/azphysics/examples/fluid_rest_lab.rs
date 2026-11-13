//! Kontrollbecken: Ruhe getrennt an beiden Wänden und in der Mitte messen.
//! cargo demo fluid_rest_lab -- --paused
#[path = "common/snapshot.rs"]
mod snapshot;

use azmath::Vector;
use azphysics::fluid::ParticleFluid;
use azphysics::geometry::Aabb;
use azviz::{Canvas, Color, Options, Plane, Scene, VizError};
use snapshot::Snapshot;
use std::path::Path;

struct RestLab {
    fluid: ParticleFluid<f32>,
    variant: u32,
    time: f64,
    steps: u64,
}

impl RestLab {
    fn new(variant: u32) -> Self {
        let full = variant == 0 || variant == 5;
        let spacing = if variant <= 1 || variant == 5 { 0.1 } else { 0.05 };
        let (width, height, fill) = if full { (1.0, 1.0, 1.0) } else { (2.0, 1.2, 0.4) };
        let mut fluid = ParticleFluid::new(
            Aabb { min: Vector([0.0, 0.0]), max: Vector([width, height]) }, spacing);
        if variant == 3 || variant == 4 {
            let layers = if variant == 3 { 1.0 } else { 3.0 };
            fluid.add_obstacle(Vector([1.0, 0.0]), Vector([1.0 + layers * spacing, 0.2]));
        }
        if variant == 5 { fluid.params.gravity = Vector::ZERO; }
        fluid.add_block(Vector([0.0, 0.0]), Vector([width, fill]));
        Self { fluid, variant, time: 0.0, steps: 0 }
    }

    fn region(&self, low: f32, high: f32) -> (usize, f32, f32) {
        let (mut count, mut sum, mut peak) = (0usize, 0.0, 0.0f32);
        for (p, v) in self.fluid.positions().iter().zip(self.fluid.velocities()) {
            if p[0] >= low && p[0] < high {
                count += 1;
                let speed = v.norm();
                sum += speed;
                peak = peak.max(speed);
            }
        }
        (count, if count == 0 { f32::NAN } else { sum / count as f32 }, peak)
    }
}

impl Scene for RestLab {
    fn title(&self) -> String { "azphysics – Ruhebecken: Diagnose".into() }
    fn update(&mut self, dt: f32) {
        self.fluid.step(dt);
        self.time += f64::from(dt);
        self.steps += 1;
    }
    fn select_variant(&mut self, variant: u32) { if variant < 6 { *self = Self::new(variant); } }
    fn state_path(&self) -> Option<&str> { Some("fluid_rest_lab.azstate") }
    fn save_state(&self, path: &Path) -> Result<(), VizError> {
        Snapshot { kind: "fluid_rest_lab_v1".into(), variant: self.variant, time: self.time, steps: self.steps,
            positions: self.fluid.positions().iter().map(|p| [p[0], p[1], 0.0]).collect(),
            velocities: self.fluid.velocities().iter().map(|p| [p[0], p[1], 0.0]).collect() }.save(path)
    }
    fn load_state(&mut self, path: &Path) -> Result<(), VizError> {
        let state = Snapshot::load(path)?;
        if state.kind != "fluid_rest_lab_v1" || state.variant >= 6
            || state.positions.iter().chain(&state.velocities).any(|p| p[2] != 0.0) {
            return Err(VizError::Encoding("Zustand gehört nicht zu diesem Versuch".into()));
        }
        let mut next = Self::new(state.variant);
        let positions: Vec<_> = state.positions.iter().map(|p| Vector([p[0], p[1]])).collect();
        let velocities: Vec<_> = state.velocities.iter().map(|p| Vector([p[0], p[1]])).collect();
        next.fluid.restore_motion(&positions, &velocities)
            .map_err(|e| VizError::Encoding(e.to_string()))?;
        next.time = state.time;
        next.steps = state.steps;
        *self = next;
        Ok(())
    }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        let bounds = self.fluid.bounds();
        let width = bounds.max[0];
        let plane = Plane::fit(canvas, [0.0, -0.1], [width, bounds.max[1] + 0.55], 24.0);
        let corners = [[0.0, bounds.max[1]], [0.0, 0.0], [width, 0.0], [width, bounds.max[1]]];
        for pair in corners.windows(2) {
            canvas.line(plane.point(pair[0]), plane.point(pair[1]), 2.0, Color::GRAY);
        }
        for obstacle in self.fluid.obstacles() { plane.rect(canvas, obstacle.min.0, obstacle.max.0, Color::GRAY); }
        for (p, v) in self.fluid.positions().iter().zip(self.fluid.velocities()) {
            let color = Color::BLUE.mix(Color::WHITE, v.norm() / 0.5);
            canvas.circle(plane.point(p.0), plane.length(self.fluid.spacing() * 0.25), color);
        }
        canvas.text([12.0, 48.0], "1 voll | 2 flach/grob | 3 flach/fein | 4 dünne Wand | 5 dicke Wand | 6 ohne g", 18.0, Color::WHITE);
        canvas.text([12.0, 72.0], &format!("Versuch {} | t={:.3} s | N={} | S speichern, L laden",
            self.variant + 1, self.time, self.fluid.positions().len()), 18.0, Color::WHITE);
        let band = 1.5 * self.fluid.spacing();
        for (row, (label, low, high)) in [("links", 0.0, band),
            ("Mitte", width * 0.35, width * 0.65), ("rechts", width - band, width + 0.001)].into_iter().enumerate() {
            let (count, mean, peak) = self.region(low, high);
            canvas.text([12.0, 96.0 + row as f32 * 22.0],
                &format!("{label}: N={count}, Mittel={mean:.6}, Maximum={peak:.6} m/s"), 18.0,
                if count > 0 && mean.is_finite() && mean < 0.02 { Color::GREEN } else { Color::ORANGE });
        }
        let finite = self.fluid.positions().iter().chain(self.fluid.velocities())
            .all(|p| p.0.iter().all(|v| v.is_finite()));
        canvas.text([12.0, 164.0], &format!("Verdichtung={:.6}, endliche Werte={finite}",
            self.fluid.mean_compression()), 18.0, if finite { Color::GRAY } else { Color::RED });
        canvas.text([12.0, 186.0], "Dichte direkt nach Laden noch nicht neu berechnet", 16.0, Color::GRAY);
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(RestLab::new(0), Options { width: 1000, height: 700, ..Options::default() })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn variants_have_nonempty_measurement_regions() {
        for variant in 0..6 {
            let scene = RestLab::new(variant);
            let width = scene.fluid.bounds().max[0];
            let band = 1.5 * scene.fluid.spacing();
            assert!(scene.region(0.0, band).0 > 0);
            assert!(scene.region(width - band, width + 0.001).0 > 0);
            assert!(scene.region(width * 0.35, width * 0.65).0 > 0);
        }
    }
}
