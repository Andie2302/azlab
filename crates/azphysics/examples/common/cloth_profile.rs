//! Gemeinsame, optionale CPU-Messung der Gewebe-Beispiele.
use azphysics::cloth::{Cloth, Collider, StepPhase};
use std::time::Instant;

#[derive(Default)]
pub struct ClothProfiler {
    pub enabled: bool,
    totals: [f64; 7],
    samples: u32,
    hits: u64,
    rollbacks: u64,
}

impl ClothProfiler {
    pub fn step(&mut self, cloth: &mut Cloth<f32>, dt: f32, colliders: &[Collider<f32>]) {
        if !self.enabled {
            cloth.step(dt, colliders);
            return;
        }
        let mut previous = Instant::now();
        let stats = cloth.step_observed(dt, colliders, &mut |phase| {
            let now = Instant::now();
            let index = match phase {
                StepPhase::Integrate => 0, StepPhase::Links => 1,
                StepPhase::Obstacles => 2, StepPhase::Proximity => 3,
                StepPhase::Repulsion => 4, StepPhase::Ccd => 5,
                StepPhase::Velocities => 6,
            };
            self.totals[index] += (now - previous).as_secs_f64() * 1000.0;
            previous = now;
        });
        self.samples += 1;
        self.hits += u64::from(stats.ccd_hits);
        self.rollbacks += u64::from(stats.full_rollbacks);
    }

    pub fn report(&mut self) -> Option<String> {
        if self.samples == 0 { return None; }
        let n = f64::from(self.samples);
        let t = self.totals.map(|x| x / n);
        let text = format!(
            "Gewebe ms/Schritt: frei {:.3}, XPBD {:.3}, Hindernisse {:.3}, Nachbarn {:.3}, Abstoßung {:.3}, CCD {:.3}, Geschwindigkeit {:.3}; CCD-Treffer {}, Rücksetzungen {} ({} Tuchschritte)",
            t[0], t[1], t[2], t[3], t[4], t[5], t[6], self.hits, self.rollbacks, self.samples);
        self.totals.fill(0.0);
        self.samples = 0;
        self.hits = 0;
        self.rollbacks = 0;
        Some(text)
    }
}
