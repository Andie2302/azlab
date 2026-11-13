//! Kontrollversuch: freier Flug ohne äußere Kräfte und Dämpfung.
//! cargo demo cloth_free_flight_no_ccd -- --paused
#[path = "common/snapshot.rs"]
mod snapshot;

use azmath::Vector;
use azphysics::cloth::Cloth;
use azviz::{Batch, Camera, Canvas, Color, Options, Scene, VizError};
use snapshot::Snapshot;
use std::path::Path;

const SIZE: usize = 12;
const VELOCITY: [f32; 3] = [0.15, 0.0, 0.05];


type Sample = [f64; 3];

/// Differenzen über vollständige Simulationsschritte, nicht über Zeichenbilder.
/// Alle Teilchen haben hier dieselbe Masse; ihr Mittel ist die Schwerpunktbewegung.
struct MotionDiagnostics {
    previous_velocity: Vec<Sample>,
    previous_acceleration: Option<Vec<Sample>>,
    acceleration: Option<(Sample, f64)>,
    jerk: Option<(Sample, f64)>,
}

impl MotionDiagnostics {
    fn new(velocities: &[Vector<f32, 3>]) -> Self {
        Self {
            previous_velocity: velocities.iter().map(|v| v.0.map(f64::from)).collect(),
            previous_acceleration: None, acceleration: None, jerk: None,
        }
    }

    fn summary(values: &[Sample]) -> (Sample, f64) {
        let mut mean = [0.0; 3];
        let mut maximum = 0.0f64;
        for value in values {
            for axis in 0..3 { mean[axis] += value[axis] / values.len() as f64; }
            maximum = maximum.max(value.iter().map(|x| x * x).sum::<f64>().sqrt());
        }
        (mean, maximum)
    }

    fn update(&mut self, velocities: &[Vector<f32, 3>], dt: f32) {
        if !dt.is_finite() || dt <= 0.0 { return; }
        let acceleration: Vec<Sample> = velocities.iter().zip(&self.previous_velocity)
            .map(|(v, old)| core::array::from_fn(|axis|
                (f64::from(v[axis]) - old[axis]) / f64::from(dt))).collect();
        self.jerk = self.previous_acceleration.as_ref().map(|old| {
            let jerk: Vec<Sample> = acceleration.iter().zip(old)
                .map(|(a, previous)| core::array::from_fn(|axis|
                    (a[axis] - previous[axis]) / f64::from(dt))).collect();
            Self::summary(&jerk)
        });
        self.acceleration = Some(Self::summary(&acceleration));
        self.previous_acceleration = Some(acceleration);
        for (old, v) in self.previous_velocity.iter_mut().zip(velocities) {
            *old = v.0.map(f64::from);
        }
    }

    fn line(name: &str, unit: &str, value: Option<(Sample, f64)>) -> String {
        match value {
            Some((mean, maximum)) => format!(
                "{name} Schwerpunkt xyz=({:+.3e}, {:+.3e}, {:+.3e}), Teilchen max={maximum:.3e} {unit}",
                mean[0], mean[1], mean[2]),
            None => format!("{name}: noch kein vollständiges Messintervall"),
        }
    }
}

struct FreeFlight {
    cloth: Cloth<f32>,
    initial: Vec<[f32; 3]>,
    time: f64,
    steps: u64,
    camera: Camera,
    diagnostics: MotionDiagnostics,
}

impl FreeFlight {
    fn new() -> Self {
        let mut cloth = Cloth::grid(SIZE, SIZE, Vector([-0.5, 0.0, -0.5]),
            Vector([1.0 / (SIZE - 1) as f32, 0.0, 0.0]),
            Vector([0.0, 0.0, 1.0 / (SIZE - 1) as f32]), 0.005);
        cloth.params.gravity = Vector::ZERO;
        cloth.params.damping = 0.0;
        cloth.params.thickness = 0.005;
        cloth.params.ccd = false; // Einziger Unterschied zum CCD-Kontrollversuch.
        cloth.set_all_velocities(Vector(VELOCITY));
        let initial = cloth.positions().iter().map(|p| p.0).collect();
        let diagnostics = MotionDiagnostics::new(cloth.velocities());
        Self { cloth, initial, diagnostics, time: 0.0, steps: 0,
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
    fn title(&self) -> String { "azphysics – Freier Flug: ohne CCD".into() }

    fn update(&mut self, dt: f32) {
        self.cloth.step(dt, &[]);
        self.diagnostics.update(self.cloth.velocities(), dt);
        self.time += f64::from(dt);
        self.steps += 1;
    }

    fn orbit(&mut self, yaw: f32, pitch: f32, zoom: f32) { self.camera.orbit(yaw, pitch, zoom); }
    fn state_path(&self) -> Option<&str> { Some("cloth_free_flight_no_ccd.azstate") }
    fn save_state(&self, path: &Path) -> Result<(), VizError> {
        Snapshot { kind: "cloth_free_flight_no_ccd_v1".into(), variant: 0, time: self.time, steps: self.steps,
            positions: self.cloth.positions().iter().map(|p| p.0).collect(),
            velocities: self.cloth.velocities().iter().map(|p| p.0).collect() }.save(path)
    }
    fn load_state(&mut self, path: &Path) -> Result<(), VizError> {
        let state = Snapshot::load(path)?;
        if state.kind != "cloth_free_flight_no_ccd_v1" || state.variant != 0 {
            return Err(VizError::Encoding("Zustand gehört nicht zu diesem Versuch".into()));
        }
        let mut next = Self::new();
        let positions: Vec<_> = state.positions.into_iter().map(Vector).collect();
        let velocities: Vec<_> = state.velocities.into_iter().map(Vector).collect();
        next.cloth.restore_motion(&positions, &velocities)
            .map_err(|e| VizError::Encoding(e.to_string()))?;
        // Ableitungen benötigen neue Intervalle; alte Historie ist nicht im Snapshot.
        next.diagnostics = MotionDiagnostics::new(next.cloth.velocities());
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
        canvas.text([12.0, 48.0], "CCD AUS | Selbstkollision AN | S speichern, L laden", 18.0, Color::WHITE);
        canvas.text([12.0, 72.0], &format!("t={:.3} s, Schritte={}, Ortsfehler={:.6} m, Geschw.-Fehler={:.6} m/s",
            self.time, self.steps, position_error, velocity_error), 18.0, color);
        canvas.text([12.0, 120.0], &MotionDiagnostics::line("a", "m/s²", self.diagnostics.acceleration), 17.0, Color::WHITE);
        canvas.text([12.0, 144.0], &MotionDiagnostics::line("Ruck", "m/s³", self.diagnostics.jerk), 17.0, Color::WHITE);
        let drift: Vec<Sample> = self.cloth.velocities().iter()
            .map(|v| core::array::from_fn(|axis| f64::from(v[axis]) - f64::from(VELOCITY[axis]))).collect();
        let (mean, _) = MotionDiagnostics::summary(&drift);
        canvas.text([12.0, 168.0], &format!(
            "Mittlerer Geschw.-Fehler xyz=({:+.3e}, {:+.3e}, {:+.3e}) m/s",
            mean[0], mean[1], mean[2]), 17.0, Color::WHITE);
        canvas.text([12.0, 192.0], "a: ab Schritt 1 | Ruck: ab Schritt 2 | nach Laden neue Messhistorie", 17.0, Color::WHITE);
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
    fn derivatives_measure_motion_and_preserve_sign() {
        let mut diagnostics = MotionDiagnostics::new(&[Vector([0.0, 0.0, 0.0]); 2]);
        assert!(diagnostics.acceleration.is_none());
        assert!(diagnostics.jerk.is_none());
        diagnostics.update(&[Vector([1.0, -2.0, 0.0]); 2], 0.5);
        assert_eq!(diagnostics.acceleration.unwrap().0, [2.0, -4.0, 0.0]);
        assert!(diagnostics.jerk.is_none());
        diagnostics.update(&[Vector([3.0, -3.0, 0.0]); 2], 0.5);
        assert_eq!(diagnostics.acceleration.unwrap().0, [4.0, -2.0, 0.0]);
        assert_eq!(diagnostics.jerk.unwrap().0, [4.0, 4.0, 0.0]);
        diagnostics.update(&[Vector([3.0, -3.0, 0.0]); 2], 0.5);
        assert_eq!(diagnostics.acceleration.unwrap().1, 0.0);
    }

    #[test]
    fn opposing_accelerations_do_not_hide_particle_maximum() {
        let mut diagnostics = MotionDiagnostics::new(&[Vector([0.0, 0.0, 0.0]); 2]);
        diagnostics.update(&[Vector([1.0, 0.0, 0.0]), Vector([-1.0, 0.0, 0.0])], 0.5);
        assert_eq!(diagnostics.acceleration.unwrap(), ([0.0; 3], 2.0));
    }

    #[test]
    fn comparison_disables_only_ccd() {
        let scene = FreeFlight::new();
        assert!(!scene.cloth.params.ccd);
        assert!(scene.cloth.params.self_collision);
        assert_eq!(scene.cloth.params.gravity, Vector::ZERO);
        assert_eq!(scene.cloth.params.damping, 0.0);
    }

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
