//! Der Zeitschritt: freie Bewegung, Abstandsbedingungen, Hindernisse.

use super::model::{Cloth, Collider, StepStats, StepPhase};
use crate::geometry::Vec3;
use azmath::Vector;
use azmath::scalar::RealElementary;

/// Verschiebt eine Lage um `amount`, ohne dass beim Runden etwas verloren
/// geht: Was die gerundete Summe nicht aufnimmt, bleibt in `carry` und wird
/// bei der nächsten Verschiebung nachgeholt (fehlerfreie Summe nach Knuth).
///
/// Lage und Rest zusammen sind die genaue Lage. Die Lage allein ist auf das
/// Raster ihres Zahltyps gerundet, und dieses Raster wird gröber, je weiter
/// das Teilchen vom Ursprung entfernt ist.
#[inline]
fn shift<T: RealElementary>(position: &mut Vec3<T>, carry: &mut Vec3<T>, amount: Vec3<T>) {
    let start = *position;
    let path = amount + *carry;
    let end = start + path;
    let taken = end - start;
    *carry = (start - (end - taken)) + (path - taken);
    *position = end;
}

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
        //
        // Der Weg eines Teilschritts ist meist nur wenige tausend
        // Rasterschritte der Lage lang; beim Addieren geht sein Rest
        // verloren, in jedem Teilschritt derselbe. Ohne Gegenmaßnahme bewegt
        // sich das Teilchen deshalb dauerhaft etwas zu schnell oder zu
        // langsam. `shift` merkt sich den Rest und holt ihn nach.
        let count = self.positions.len();
        self.carry.resize(count, Vector::ZERO);
        self.predicted.resize(count, Vector::ZERO);
        self.predicted_carry.resize(count, Vector::ZERO);
        let keep = (T::ONE - self.params.damping * h).max(T::ZERO);
        for index in 0..count {
            self.previous[index] = self.positions[index];
            if self.inverse_mass[index] > T::ZERO {
                let velocity = (self.velocities[index] + self.params.gravity * h) * keep;
                self.velocities[index] = velocity;
                shift(&mut self.positions[index], &mut self.carry[index], velocity * h);
            }
            self.predicted[index] = self.positions[index];
            self.predicted_carry[index] = self.carry[index];
        }

        observe(StepPhase::Integrate);

        // 2. Abstandsbedingungen (XPBD).
        self.lambdas.fill(T::ZERO);
        // Die Durchläufe gehen abwechselnd vorwärts und rückwärts durch die
        // Verbindungen, immer paarweise. Läuft man stets in derselben
        // Richtung, ist die Korrektur als Ganzes unsymmetrisch: Sie dreht
        // Fehler im Kreis, statt sie nur zu verkleinern, und zusammen mit
        // der Geschwindigkeit aus der Lageänderung schaukelt sich das auf.
        // Ein ruhig hängendes Tuch gerät so von selbst in Bewegung. Hin und
        // zurück zusammen sind symmetrisch und verkleinern jeden Fehler.
        let sweeps = self.params.iterations.div_ceil(2) * 2;
        for sweep in 0..sweeps {
            self.solve_links(h, sweep % 2 == 1);
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

        // 4. Geschwindigkeit: die der freien Bewegung, geändert um das, was
        // die Korrekturen des Teilschritts an der Lage verschoben haben.
        //
        // Die Geschwindigkeit wird nicht aus dem ganzen zurückgelegten Weg
        // zurückgerechnet. Der ist auf das Raster der Lage gerundet, und die
        // Geschwindigkeit würde in jedem Teilschritt auf dieses Raster
        // einrasten: ein Fehler, der mit dem Abstand vom Ursprung wächst und
        // ein frei fliegendes Tuch stetig abtreiben lässt. Ein Teilchen, das
        // nichts korrigiert hat, behält seine Geschwindigkeit so bitgenau.
        for index in 0..self.positions.len() {
            if self.inverse_mass[index] <= T::ZERO {
                continue;
            }
            let mut velocity = if self.positions[index] == self.previous[index] {
                // An den Anfang des Teilschritts zurückgesetzt: steht.
                Vector::ZERO
            } else {
                let moved = (self.positions[index] - self.predicted[index])
                    + (self.carry[index] - self.predicted_carry[index]);
                self.velocities[index] + moved / h
            };
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
    fn solve_links(&mut self, h: T, backward: bool) {
        let inverse_h2 = (h * h).recip();
        let count = self.links.len();
        for step in 0..count {
            let index = if backward { count - 1 - step } else { step };
            let (link, lambda) = (&self.links[index], &mut self.lambdas[index]);
            let (wi, wj) = (self.inverse_mass[link.i], self.inverse_mass[link.j]);
            if wi + wj <= T::ZERO {
                continue;
            }
            // Abstand aus den genauen Lagen (Lage plus Rest). Aus den
            // gerundeten Lagen allein erschiene ein Tuch, das als Ganzes
            // fliegt, um einen Rasterschritt verzerrt, und die Korrektur
            // dieses Scheinfehlers verteilte Stöße in Rastergröße.
            let delta = (self.positions[link.i] - self.positions[link.j]) + (self.carry[link.i] - self.carry[link.j]);
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
            shift(&mut self.positions[link.i], &mut self.carry[link.i], direction * (change * wi));
            shift(&mut self.positions[link.j], &mut self.carry[link.j], -direction * (change * wj));
        }
    }

    /// Schiebt Teilchen aus einem Hindernis und bremst ihr Gleiten.
    fn collide(&mut self, collider: &Collider<T>) {
        let thickness = self.params.thickness;
        let friction = self.params.friction.max(T::ZERO);
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
            // Reibung nach Coulomb. Die Eindringtiefe steht für die Kraft,
            // mit der das Teilchen auf die Oberfläche drückt; das Gleiten
            // längs der Oberfläche darf um höchstens Reibwert mal Tiefe
            // gekürzt werden. Gleitet das Teilchen weniger weit, haftet es.
            // Ruhend auf einer Schräge heißt das: Es bleibt liegen, solange
            // die Steigung kleiner ist als der Reibwert – unabhängig von der
            // Schrittweite.
            let movement = lifted - self.previous[index];
            let along = movement - normal * movement.dot(normal);
            let slip = along.norm();
            let grip = friction * depth;
            self.positions[index] =
                if slip <= grip || slip <= T::ZERO { lifted - along } else { lifted - along * (grip / slip) };
        }
    }
}
