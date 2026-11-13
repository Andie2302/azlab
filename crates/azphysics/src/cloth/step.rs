//! Der Zeitschritt: freie Bewegung, Abstandsbedingungen, Hindernisse.

use super::model::{Cloth, Collider, StepStats, StepPhase};
use azmath::Vector;
use azmath::scalar::RealElementary;

impl<T: RealElementary> Cloth<T> {
    /// Rechnet die Simulation um `dt` Sekunden weiter.
    ///
    /// Der Schritt wird in `params.substeps` Teilschritte zerlegt. Für eine
    /// Bildfolge sind `dt = 1/60` und acht bis zwanzig Teilschritte üblich.
    pub fn step(&mut self, dt: T, colliders: &[Collider<T>]) -> StepStats {
        self.step_observed(dt, colliders, &mut |_| {})
    }

    /// Wie `step`, ruft nach jedem Abschnitt `observe` auf.
    /// Der Callback darf messen, verändert aber weder Zeitschritt noch Physik.
    /// Ohne Messbedarf eliminiert der Compiler den leeren Callback.
    pub fn step_observed(
        &mut self, dt: T, colliders: &[Collider<T>],
        observe: &mut impl FnMut(StepPhase),
    ) -> StepStats {
        let substeps = self.params.substeps.max(1);
        let h = dt / T::from_f64(f64::from(substeps));
        let mut stats = StepStats::default();
        if !(h > T::ZERO && h.is_finite()) {
            return stats;
        }
        for _ in 0..substeps {
            self.substep(h, colliders, &mut stats, observe);
        }
        stats
    }

    fn substep(&mut self, h: T, colliders: &[Collider<T>], stats: &mut StepStats, observe: &mut impl FnMut(StepPhase)) {
        // 1. Freie Bewegung unter Schwerkraft.
        let keep = (T::ONE - self.params.damping * h).max(T::ZERO);
        for index in 0..self.positions.len() {
            self.previous[index] = self.positions[index];
            if self.inverse_mass[index] > T::ZERO {
                let velocity = (self.velocities[index] + self.params.gravity * h) * keep;
                self.velocities[index] = velocity;
                self.positions[index] += velocity * h;
            }
        }

        observe(StepPhase::Integrate);

        // 2. Abstandsbedingungen (XPBD).
        self.lambdas.fill(T::ZERO);
        for _ in 0..self.params.iterations.max(1) {
            self.solve_links(h);
        }

        observe(StepPhase::Links);

        // 3. Hindernisse.
        for collider in colliders {
            self.collide(collider);
        }
        observe(StepPhase::Obstacles);
        self.pushes.clear();
        self.pushes.resize(self.positions.len(), Vector::ZERO);
        if self.params.self_collision {
            self.repel(stats, observe);
            observe(StepPhase::Repulsion);
        }
        if self.params.ccd {
            self.prevent_penetration(stats, observe);
            observe(StepPhase::Ccd);
        }

        // 4. Geschwindigkeit aus der tatsächlich erreichten Lage.
        for index in 0..self.positions.len() {
            if self.inverse_mass[index] <= T::ZERO {
                continue;
            }
            let mut velocity = (self.positions[index] - self.previous[index]) / h;
            // Abstoßung ist unelastisch: Sie darf die Annäherung aufhalten,
            // aber das Teilchen nicht wegschleudern. Was die Verschiebung an
            // Geschwindigkeit vom Dreieck weg hinzugefügt hat, entfällt.
            let push = self.pushes[index];
            if let Some(away) = push.normalize() {
                let gained = velocity.dot(away);
                let allowed = (gained - push.norm() / h).max(T::ZERO);
                if gained > allowed {
                    velocity -= away * (gained - allowed);
                }
            }
            self.velocities[index] = velocity;
        }
        observe(StepPhase::Velocities);
    }

    /// Ein Durchlauf über alle Abstandsbedingungen.
    fn solve_links(&mut self, h: T) {
        let inverse_h2 = (h * h).recip();
        for (link, lambda) in self.links.iter().zip(self.lambdas.iter_mut()) {
            let (wi, wj) = (self.inverse_mass[link.i], self.inverse_mass[link.j]);
            if wi + wj <= T::ZERO {
                continue;
            }
            let delta = self.positions[link.i] - self.positions[link.j];
            let length = delta.norm();
            if length <= T::ZERO || length.is_nan() {
                continue;
            }
            let compliance = if link.bend {
                self.params.bend_compliance
            } else {
                self.params.stretch_compliance
            };
            let alpha = compliance * inverse_h2;
            let change = (link.rest - length - alpha * *lambda) / (wi + wj + alpha);
            *lambda += change;
            let direction = delta / length;
            self.positions[link.i] += direction * (change * wi);
            self.positions[link.j] -= direction * (change * wj);
        }
    }

    /// Schiebt Teilchen aus einem Hindernis und bremst ihr Gleiten.
    fn collide(&mut self, collider: &Collider<T>) {
        let thickness = self.params.thickness;
        let slide = (T::ONE - self.params.friction).max(T::ZERO).min(T::ONE);
        for index in 0..self.positions.len() {
            if self.inverse_mass[index] <= T::ZERO {
                continue;
            }
            let position = self.positions[index];
            // Normale am Berührpunkt und Eindringtiefe (positiv = eingedrungen).
            let contact = match *collider {
                Collider::Plane { point, normal } => {
                    let depth = thickness - (position - point).dot(normal);
                    (depth > T::ZERO).then_some((normal, depth))
                }
                Collider::Sphere { center, radius } => {
                    let offset = position - center;
                    let distance = offset.norm();
                    let depth = radius + thickness - distance;
                    if depth > T::ZERO {
                        let normal = if distance > T::ZERO {
                            offset / distance
                        } else {
                            Vector([T::ZERO, T::ONE, T::ZERO])
                        };
                        Some((normal, depth))
                    } else {
                        None
                    }
                }
            };
            let Some((normal, depth)) = contact else {
                continue;
            };
            let lifted = position + normal * depth;
            // Reibung: den Anteil der Bewegung längs der Oberfläche kürzen.
            let movement = lifted - self.previous[index];
            let along = movement - normal * movement.dot(normal);
            self.positions[index] = lifted - along * (T::ONE - slide);
        }
    }
}
