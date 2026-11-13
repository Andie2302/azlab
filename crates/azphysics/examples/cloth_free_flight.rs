//! Kontrollversuch: freier Flug ohne äußere Kräfte und Dämpfung.
//! cargo demo cloth_free_flight -- --paused
#[path = "common/snapshot.rs"]
mod snapshot;

use azmath::Vector;
use azphysics::cloth::Cloth;
use azviz::{Batch, Camera, Canvas, Color, Options, Scene, VizError};
use snapshot::Snapshot;
use std::path::Path;

const SIZE: usize = 12;
const VELOCITY: [f32; 3] = [0.15, 0.0, 0.05];

struct FreeFlight {
    cloth: Cloth<f32>,
    initial: Vec<[f32; 3]>,
    time: f64,
    steps: u64,
    camera: Camera,
}

impl FreeFlight {
    fn new() -> Self {
        let mut cloth = Cloth::grid(SIZE, SIZE, Vector([-0.5, 0.0, -0.5]),
            Vector([1.0 / (SIZE - 1) as f32, 0.0, 0.0]),
            Vector([0.0, 0.0, 1.0 / (SIZE - 1) as f32]), 0.005);
        cloth.params.gravity = Vector::ZERO;
        cloth.params.damping = 0.0;
        cloth.params.thickness = 0.005;
        cloth.set_all_velocities(Vector(VELOCITY));
        let initial = cloth.positions().iter().map(|p| p.0).collect();
        Self { cloth, initial, time: 0.0, steps: 0,
            camera: Camera::look_at([1.5, 1.2, 2.0], [0.0, 0.0, 0.0], 40.0) }
    }

    fn errors(&self) -> (f32, f32) {
        let offset = Vector(VELOCITY) * self.time as f32;
        let position = self.cloth.positions().iter().zip(&self.initial)
            .map(|(p, &start)| (*p - Vector(start) - offset).norm()).fold(0.0f32, f32::max);
        let velocity = self.cloth.velocities().iter()
            .map(|v| (*v - Vector(VELOCITY)).norm()).fold(0.0f32, f32::max);
        (position, velocity)
    }
}

impl Scene for FreeFlight {
    fn title(&self) -> String { "azphysics – Freier Flug: Kontrollversuch".into() }

    fn update(&mut self, dt: f32) {
        self.cloth.step(dt, &[]);
        self.time += f64::from(dt);
        self.steps += 1;
    }

    fn orbit(&mut self, yaw: f32, pitch: f32, zoom: f32) { self.camera.orbit(yaw, pitch, zoom); }
    fn state_path(&self) -> Option<&str> { Some("cloth_free_flight.azstate") }
    fn save_state(&self, path: &Path) -> Result<(), VizError> {
        Snapshot { kind: "cloth_free_flight_v1".into(), variant: 0, time: self.time, steps: self.steps,
            positions: self.cloth.positions().iter().map(|p| p.0).collect(),
            velocities: self.cloth.velocities().iter().map(|p| p.0).collect() }.save(path)
    }
    fn load_state(&mut self, path: &Path) -> Result<(), VizError> {
        let state = Snapshot::load(path)?;
        if state.kind != "cloth_free_flight_v1" || state.variant != 0 {
            return Err(VizError::Encoding("Zustand gehört nicht zu diesem Versuch".into()));
        }
        let mut next = Self::new();
        let positions: Vec<_> = state.positions.into_iter().map(Vector).collect();
        let velocities: Vec<_> = state.velocities.into_iter().map(Vector).collect();
        next.cloth.restore_motion(&positions, &velocities)
            .map_err(|e| VizError::Encoding(e.to_string()))?;
        next.time = state.time;
        next.steps = state.steps;
        next.camera = self.camera;
        *self = next;
        Ok(())
    }

    fn draw(&self, canvas: &mut dyn Canvas) {
        canvas.clear(Color::BACKGROUND);
        // Mitbewegte Darstellung: gleichförmiger Flug bleibt im Bild.
        let offset = Vector(VELOCITY) * self.time as f32;
        let positions: Vec<_> = self.cloth.positions().iter().map(|p| (*p - offset).0).collect();
        let (position_error, velocity_error) = self.errors();
        let finite = self.cloth.positions().iter().chain(self.cloth.velocities())
            .all(|p| p.0.iter().all(|x| x.is_finite()));
        let color = if finite && position_error < 0.001 && velocity_error < 0.001 { Color::GREEN } else { Color::RED };
        let mut batch = Batch::new();
        batch.mesh(&positions, self.cloth.triangles(), color, color);
        batch.draw(canvas, &self.camera, true);
        canvas.text([12.0, 48.0], "Mitbewegte Ansicht | keine Schwerkraft, keine Dämpfung | S speichern, L laden", 18.0, Color::WHITE);
        canvas.text([12.0, 72.0], &format!("t={:.3} s, Schritte={}, Ortsfehler={:.6} m, Geschw.-Fehler={:.6} m/s",
            self.time, self.steps, position_error, velocity_error), 18.0, color);
        canvas.text([12.0, 96.0], &format!("Dehnung={:.6}, endliche Werte={finite}", self.cloth.max_strain()), 18.0, color);
    }
}

fn main() -> Result<(), VizError> {
    azviz::run(FreeFlight::new(), Options { width: 1000, ..Options::default() })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unforced_motion_matches_translation() {
        let mut scene = FreeFlight::new();
        for _ in 0..120 { scene.update(1.0 / 60.0); }
        let (p, v) = scene.errors();
        assert!(p < 0.001, "Ortsfehler {p}");
        assert!(v < 0.001, "Geschwindigkeitsfehler {v}");
        assert!(scene.cloth.max_strain() < 0.001);
    }
}
