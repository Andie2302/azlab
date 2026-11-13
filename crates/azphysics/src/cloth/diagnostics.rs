//! Messwerte aus den öffentlich sichtbaren (gerundeten) Teilchenzuständen.

use super::Cloth;
use azmath::scalar::RealElementary;

/// Momentaufnahme des beweglichen Tuchs, ohne Speicherallokation.
///
/// Festgehaltene Teilchen und Hindernisse haben inverse Masse null und gehen
/// nicht in Masse, Energie oder Impuls ein. Haltepunkte können daher äußere
/// Kräfte auf dieses Teilsystem ausüben. Summiert wird in `f64`, auch bei
/// `Cloth<f32>`; Wurzeln verwendet der Rechen-Backend des Tuchs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClothDiagnostics {
    /// Alle Eingaben und berechneten Messwerte sind endlich und gültig.
    pub finite: bool,
    /// Anzahl beweglicher Teilchen.
    pub moving_particles: usize,
    /// Summe der beweglichen Massen in kg.
    pub mass: f64,
    /// Schwerpunkt in m; bei Masse null der Nullvektor.
    pub center_of_mass: [f64; 3],
    /// Kinetische Energie in J.
    pub kinetic_energy: f64,
    /// Gravitationsenergie `-Σ m g·x` in J, Bezugspunkt Weltursprung.
    pub gravitational_energy: f64,
    /// Energie nachgiebiger Verbindungen `Σ C²/(2 compliance)` in J.
    /// Harte Bedingungen (compliance null) tragen keine Federenergie bei.
    /// Bei abgeschalteten Verbindungen (`iterations = 0`) ist sie null.
    pub elastic_energy: f64,
    /// Summe aus kinetischer, Gravitations- und Federenergie in J.
    /// Kontaktenergie ist nicht enthalten; Reibung, Dämpfung und unelastische
    /// Kollisionen erlauben Energieverlust.
    pub mechanical_energy: f64,
    /// Linearer Impuls in kg m/s, Komponenten x/y/z.
    pub momentum: [f64; 3],
    /// Drehimpuls um den momentanen Schwerpunkt in kg m²/s.
    pub angular_momentum: [f64; 3],
    /// Größte absolute relative Längenabweichung der Dehn-/Scherbedingungen.
    pub max_strain: f64,
    /// RMS der relativen Längenabweichungen der Dehn-/Scherbedingungen.
    pub rms_strain: f64,
    /// Größte absolute Abstandsverletzung aller Verbindungen in m,
    /// einschließlich der als Abstandsbedingungen modellierten Biegung.
    pub max_constraint_error: f64,
}

impl<T: RealElementary> Cloth<T> {
    /// Misst den aktuellen Zustand, ohne ihn oder die Simulation zu ändern.
    ///
    /// `finite` muss vor der Interpretation geprüft werden: Insbesondere
    /// dürfen NaN-Werte nicht durch Maximumsbildung unbemerkt verschwinden.
    /// Energie- und Impulserhaltung sind nur bei passenden Randbedingungen
    /// zu erwarten, nicht bei Haltepunkten, Dämpfung oder äußeren Kontakten.
    pub fn diagnostics(&self) -> ClothDiagnostics {
        let mut result = ClothDiagnostics {
            finite: true,
            moving_particles: 0,
            mass: 0.0,
            center_of_mass: [0.0; 3],
            kinetic_energy: 0.0,
            gravitational_energy: 0.0,
            elastic_energy: 0.0,
            mechanical_energy: 0.0,
            momentum: [0.0; 3],
            angular_momentum: [0.0; 3],
            max_strain: 0.0,
            rms_strain: 0.0,
            max_constraint_error: 0.0,
        };
        let gravity = self.params.gravity.0.map(|x| x.to_f64());
        result.finite &= gravity.iter().all(|x| x.is_finite());
        for ((position, velocity), inverse) in self
            .positions
            .iter()
            .zip(&self.velocities)
            .zip(&self.inverse_mass)
        {
            let (p, v, w) = (
                position.0.map(|x| x.to_f64()),
                velocity.0.map(|x| x.to_f64()),
                inverse.to_f64(),
            );
            result.finite &= p.iter().chain(&v).all(|x| x.is_finite()) && w.is_finite() && w >= 0.0;
            if w <= 0.0 {
                continue;
            }
            let mass = 1.0 / w;
            result.mass += mass;
            result.moving_particles += 1;
            for axis in 0..3 {
                result.center_of_mass[axis] += mass * p[axis];
                result.momentum[axis] += mass * v[axis];
                result.kinetic_energy += 0.5 * mass * v[axis] * v[axis];
                result.gravitational_energy -= mass * gravity[axis] * p[axis];
            }
        }
        if result.mass > 0.0 {
            for x in &mut result.center_of_mass {
                *x /= result.mass;
            }
        }
        for ((position, velocity), inverse) in self
            .positions
            .iter()
            .zip(&self.velocities)
            .zip(&self.inverse_mass)
        {
            let w = inverse.to_f64();
            if w <= 0.0 {
                continue;
            }
            let r: [f64; 3] =
                core::array::from_fn(|i| position[i].to_f64() - result.center_of_mass[i]);
            let p = velocity.0.map(|x| x.to_f64() / w);
            for (axis, a, b) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
                result.angular_momentum[axis] += r[a] * p[b] - r[b] * p[a];
            }
        }
        let mut strain_squares = 0.0;
        let mut strain_count = 0;
        for link in &self.links {
            let rest = link.rest.to_f64();
            // Abstände zuerst in f64 bilden: gemeinsame Translation darf
            // nicht zusätzlich durch eine f32-Subtraktion verfälscht werden.
            let squared: f64 = (0..3)
                .map(|axis| {
                    let d = self.positions[link.i][axis].to_f64()
                        - self.positions[link.j][axis].to_f64();
                    d * d
                })
                .sum();
            let length = T::from_f64(squared).sqrt().to_f64();
            let error = length - rest;
            let compliance = if link.bend {
                self.params.bend_compliance
            } else {
                self.params.stretch_compliance
            }
            .to_f64();
            result.finite &=
                error.is_finite() && rest >= 0.0 && compliance.is_finite() && compliance >= 0.0;
            result.max_constraint_error = result.max_constraint_error.max(error.abs());
            if !link.bend && rest > 0.0 {
                let strain = error / rest;
                result.max_strain = result.max_strain.max(strain.abs());
                strain_squares += strain * strain;
                strain_count += 1;
            }
            if self.params.iterations > 0 && compliance > 0.0 {
                result.elastic_energy += error * error / (2.0 * compliance);
            }
        }
        if strain_count > 0 {
            result.rms_strain = T::from_f64(strain_squares / f64::from(strain_count))
                .sqrt()
                .to_f64();
        }
        result.mechanical_energy =
            result.kinetic_energy + result.gravitational_energy + result.elastic_energy;
        result.finite &= [
            result.mass,
            result.kinetic_energy,
            result.gravitational_energy,
            result.elastic_energy,
            result.mechanical_energy,
            result.max_strain,
            result.rms_strain,
            result.max_constraint_error,
        ]
        .iter()
        .chain(&result.center_of_mass)
        .chain(&result.momentum)
        .chain(&result.angular_momentum)
        .all(|x| x.is_finite());
        result
    }
}

#[cfg(all(test, any(feature = "std", feature = "libm")))]
mod tests {
    use super::*;
    use azmath::Vector;

    fn cloth() -> Cloth<f64> {
        let mut c = Cloth::grid(
            2,
            2,
            Vector([0.0; 3]),
            Vector([1.0, 0.0, 0.0]),
            Vector([0.0, 0.0, 1.0]),
            2.0,
        );
        c.params.gravity = Vector([0.0, -10.0, 0.0]);
        c
    }

    #[test]
    fn translation_mass_energy_and_obstacles() {
        let mut c = cloth();
        c.set_all_velocities(Vector([3.0, -2.0, 1.0]));
        c.add_obstacle_triangle(Vector([100.0; 3]), Vector([101.0; 3]), Vector([102.0; 3]));
        let d = c.diagnostics();
        assert!(d.finite);
        assert_eq!(d.mass, 8.0);
        assert_eq!(d.moving_particles, 4);
        assert_eq!(d.center_of_mass, [0.5, 0.0, 0.5]);
        assert_eq!(d.momentum, [24.0, -16.0, 8.0]);
        assert_eq!(d.kinetic_energy, 56.0);
        assert_eq!(d.angular_momentum, [0.0; 3]);
        c.pin(0);
        assert_eq!(c.diagnostics().mass, 6.0);
    }

    #[test]
    fn rotation_is_independent_of_translation() {
        let mut c = cloth();
        for i in 0..4 {
            let p = c.positions()[i];
            c.set_velocity(i, Vector([p[2] - 0.5, 0.0, -(p[0] - 0.5)]));
        }
        let d = c.diagnostics();
        assert_eq!(d.momentum, [0.0; 3]);
        assert_eq!(d.angular_momentum, [0.0, 4.0, 0.0]);
        for i in 0..4 {
            c.set_position(i, c.positions()[i] + Vector([10.0, 3.0, -7.0]));
            c.set_velocity(i, c.velocities()[i] + Vector([2.0, 5.0, -3.0]));
        }
        let shifted = c.diagnostics();
        assert_eq!(shifted.angular_momentum, d.angular_momentum);
        assert_eq!(shifted.gravitational_energy, 240.0);
    }

    #[test]
    fn constraint_energy_and_invalid_values() {
        let mut c = cloth();
        c.params.stretch_compliance = 0.1;
        c.set_position(0, Vector([-0.25, 0.0, 0.0]));
        let d = c.diagnostics();
        assert!(d.finite && d.elastic_energy > 0.0 && d.rms_strain > 0.0);
        assert!((d.max_strain - 0.25).abs() < 1e-12);
        assert!((d.max_constraint_error - 0.25).abs() < 1e-12);
        c.set_position(0, Vector([f64::NAN, 0.0, 0.0]));
        assert!(!c.diagnostics().finite);
    }

    #[test]
    fn fixed_system_has_no_moving_mass() {
        let mut c = cloth();
        for i in 0..4 {
            c.pin(i);
        }
        let d = c.diagnostics();
        assert!(d.finite);
        assert_eq!(d.mass, 0.0);
        assert_eq!(d.center_of_mass, [0.0; 3]);
        assert_eq!(d.mechanical_energy, 0.0);
    }
}
