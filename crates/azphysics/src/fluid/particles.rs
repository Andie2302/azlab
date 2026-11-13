//! Flüssigkeit aus Teilchen (*position based fluids*).

use super::Vec2;
use crate::geometry::Aabb;
use alloc::vec::Vec;
use azmath::Vector;
use azmath::scalar::RealElementary;

/// Einstellungen der Teilchenflüssigkeit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParticleParams<T> {
    /// Schwerebeschleunigung.
    pub gravity: Vec2<T>,
    /// Teilschritte je Zeitschritt.
    pub substeps: u32,
    /// Durchläufe der Dichtekorrektur je Teilschritt.
    pub iterations: u32,
    /// Weichheit der Dichtebedingung; größer ist nachgiebiger und ruhiger.
    pub relaxation: T,
    /// Zähigkeit: gleicht die Geschwindigkeit benachbarter Teilchen an
    /// (0 bis 1).
    pub viscosity: T,
    /// Stärke der Abstoßung auf kurze Entfernung; verhindert, dass Teilchen
    /// verklumpen.
    pub anti_clumping: T,
}

impl<T: RealElementary> Default for ParticleParams<T> {
    fn default() -> Self {
        Self {
            gravity: Vector([T::ZERO, T::from_f64(-9.81)]),
            substeps: 2,
            iterations: 4,
            relaxation: T::from_f64(1.0),
            viscosity: T::from_f64(0.05),
            anti_clumping: T::from_f64(0.002),
        }
    }
}

/// Reichweite der Kerne in Teilchenabständen.
const SUPPORT: f64 = 2.4;

/// Dicke der Wände in Teilchenlagen.
const WALL_LAYERS: f64 = 3.0;

/// Flüssigkeit mit freier Oberfläche in einem rechteckigen Becken.
///
/// Verfahren nach Macklin und Müller, „Position Based Fluids“ (2013): Die
/// Dichte an jedem Teilchen ergibt sich aus den Nachbarn im Umkreis; liegt
/// sie über der Ruhedichte, werden die Teilchen auseinandergeschoben. Die
/// Wände und Hindernisse bestehen aus ruhenden Teilchen, die zur Dichte
/// beitragen und selbst Druck aufbauen, wenn die Flüssigkeit sich an ihnen
/// staut.
#[derive(Clone, Debug)]
pub struct ParticleFluid<T> {
    bounds: Aabb<T, 2>,
    spacing: T,
    /// Reichweite der Kerne.
    support: T,
    rest_density: T,
    positions: Vec<Vec2<T>>,
    previous: Vec<Vec2<T>>,
    velocities: Vec<Vec2<T>>,
    walls: Vec<Vec2<T>>,
    obstacles: Vec<Aabb<T, 2>>,
    densities: Vec<T>,
    // Arbeitsspeicher.
    lambdas: Vec<T>,
    deltas: Vec<Vec2<T>>,
    cell_start: Vec<usize>,
    cell_items: Vec<usize>,
    neighbour_start: Vec<usize>,
    neighbours: Vec<usize>,
    /// Einstellungen; dürfen zwischen den Schritten geändert werden.
    pub params: ParticleParams<T>,
}

impl<T: RealElementary> ParticleFluid<T> {
    /// Leeres Becken `bounds`; die Teilchen haben den Abstand `spacing`.
    pub fn new(bounds: Aabb<T, 2>, spacing: T) -> Self {
        let spacing = if spacing > T::ZERO && spacing.is_finite() { spacing } else { T::ONE };
        let support = spacing * T::from_f64(SUPPORT);
        let mut fluid = Self {
            bounds,
            spacing,
            support,
            rest_density: T::ONE,
            positions: Vec::new(),
            previous: Vec::new(),
            velocities: Vec::new(),
            walls: Vec::new(),
            obstacles: Vec::new(),
            densities: Vec::new(),
            lambdas: Vec::new(),
            deltas: Vec::new(),
            cell_start: Vec::new(),
            cell_items: Vec::new(),
            neighbour_start: Vec::new(),
            neighbours: Vec::new(),
            params: ParticleParams::default(),
        };
        // Ruhedichte: Dichte in einem vollbesetzten quadratischen Gitter.
        let reach = SUPPORT as i32 + 1;
        let mut density = T::ZERO;
        for a in -reach..=reach {
            for b in -reach..=reach {
                let offset: Vec2<T> = Vector([T::from_f64(f64::from(a)), T::from_f64(f64::from(b))]) * spacing;
                density += fluid.density_kernel(offset.norm());
            }
        }
        fluid.rest_density = density;
        // Wände: drei Lagen ruhender Teilchen rund um das Becken. So viele
        // braucht die innerste Lage, um selbst die volle Dichte zu sehen.
        let outer = bounds.expanded(spacing * T::from_f64(WALL_LAYERS));
        fluid.lattice(outer.min, outer.max, |point, fluid| {
            if !bounds.contains(point) {
                fluid.walls.push(point);
            }
        });
        fluid
    }

    /// Ruft `visit` für jeden Punkt des Teilchengitters im Rechteck auf.
    fn lattice(&mut self, min: Vec2<T>, max: Vec2<T>, mut visit: impl FnMut(Vec2<T>, &mut Self)) {
        let half = self.spacing * T::from_f64(0.5);
        // Gitter am Becken ausgerichtet, damit Wände, Hindernisse und
        // Flüssigkeit lückenlos aneinanderpassen.
        let snap = |value: T, origin: T| ((value - origin) / self.spacing).floor().to_f64();
        let first = [snap(min[0], self.bounds.min[0]), snap(min[1], self.bounds.min[1])];
        let last = [snap(max[0] - half, self.bounds.min[0]), snap(max[1] - half, self.bounds.min[1])];
        if !(first.iter().chain(&last).all(|value| value.is_finite())) {
            return;
        }
        let mut b = first[1];
        while b <= last[1] {
            let mut a = first[0];
            while a <= last[0] {
                let point = Vector([
                    self.bounds.min[0] + T::from_f64(a) * self.spacing + half,
                    self.bounds.min[1] + T::from_f64(b) * self.spacing + half,
                ]);
                if point[0] >= min[0] && point[0] <= max[0] && point[1] >= min[1] && point[1] <= max[1] {
                    visit(point, self);
                }
                a += 1.0;
            }
            b += 1.0;
        }
    }

    /// Füllt das Rechteck mit ruhender Flüssigkeit; liefert die Zahl der
    /// neuen Teilchen. Was außerhalb des Beckens oder in einem Hindernis
    /// läge, entfällt.
    pub fn add_block(&mut self, min: Vec2<T>, max: Vec2<T>) -> usize {
        let before = self.positions.len();
        self.lattice(min, max, |point, fluid| {
            if fluid.bounds.contains(point) && !fluid.obstacles.iter().any(|block| block.contains(point)) {
                fluid.positions.push(point);
            }
        });
        let count = self.positions.len();
        self.previous.resize(count, Vector::ZERO);
        self.velocities.resize(count, Vector::ZERO);
        self.densities.resize(count, T::ONE);
        self.deltas.resize(count, Vector::ZERO);
        count - before
    }

    /// Stellt ein festes rechteckiges Hindernis in das Becken.
    ///
    /// Es sollte mindestens drei Teilchenabstände dick sein: Dünnere
    /// Hindernisse bauen nicht den vollen Gegendruck auf, und das Wasser
    /// kommt an ihnen nicht ganz zur Ruhe.
    pub fn add_obstacle(&mut self, min: Vec2<T>, max: Vec2<T>) {
        self.lattice(min, max, |point, fluid| fluid.walls.push(point));
        self.obstacles.push(Aabb { min, max });
    }

    /// Das Becken.
    #[inline]
    pub fn bounds(&self) -> Aabb<T, 2> {
        self.bounds
    }

    /// Abstand der Teilchen in Ruhe.
    #[inline]
    pub fn spacing(&self) -> T {
        self.spacing
    }

    /// Lagen der Flüssigkeitsteilchen.
    #[inline]
    pub fn positions(&self) -> &[Vec2<T>] {
        &self.positions
    }

    /// Geschwindigkeiten der Flüssigkeitsteilchen.
    #[inline]
    pub fn velocities(&self) -> &[Vec2<T>] {
        &self.velocities
    }

    /// Stellt Bewegung im bestehenden Becken wieder her. Geometrie und
    /// Parameter müssen zur Aufzeichnung passen. Dichten werden im nächsten
    /// Schritt neu berechnet; direkt nach dem Laden stehen sie auf eins.
    /// Bei ungültigen Daten bleibt der Zustand vollständig erhalten.
    pub fn restore_motion(&mut self, positions: &[Vec2<T>], velocities: &[Vec2<T>])
        -> Result<(), crate::MotionError>
    {
        if positions.len() != self.positions.len() || velocities.len() != self.velocities.len() {
            return Err(crate::MotionError::DimensionMismatch);
        }
        if positions.iter().chain(velocities).any(|v| v.0.iter().any(|x| !x.is_finite())) {
            return Err(crate::MotionError::NonFinite);
        }
        if positions.iter().any(|&p| !self.bounds.contains(p)
            || self.obstacles.iter().any(|b| b.contains(p))) {
            return Err(crate::MotionError::InvalidPosition);
        }
        self.positions.copy_from_slice(positions);
        self.previous.copy_from_slice(positions);
        self.velocities.copy_from_slice(velocities);
        self.densities.fill(T::ONE);
        self.neighbour_start.clear();
        self.neighbours.clear();
        Ok(())
    }

    /// Dichte an jedem Teilchen, bezogen auf die Ruhedichte.
    #[inline]
    pub fn densities(&self) -> &[T] {
        &self.densities
    }

    /// Die ruhenden Teilchen der Wände und Hindernisse.
    #[inline]
    pub fn walls(&self) -> &[Vec2<T>] {
        &self.walls
    }

    /// Die Hindernisse im Becken.
    #[inline]
    pub fn obstacles(&self) -> &[Aabb<T, 2>] {
        &self.obstacles
    }

    /// Mittlere Verdichtung über die Ruhedichte hinaus; null heißt
    /// nirgends zusammengedrückt.
    pub fn mean_compression(&self) -> T {
        if self.densities.is_empty() {
            return T::ZERO;
        }
        let sum = self.densities.iter().fold(T::ZERO, |sum, &density| sum + (density - T::ONE).max(T::ZERO));
        sum / T::from_f64(self.densities.len() as f64)
    }

    /// Kern für die Dichte (poly6, in der Ebene normiert).
    #[inline]
    fn density_kernel(&self, distance: T) -> T {
        let q = distance / self.support;
        if q >= T::ONE {
            return T::ZERO;
        }
        let base = T::ONE - q * q;
        T::from_f64(4.0) / (T::PI * self.support * self.support) * base * base * base
    }

    /// Betrag der Ableitung des Kerns für die Korrektur (spiky).
    #[inline]
    fn gradient_kernel(&self, distance: T) -> T {
        let q = distance / self.support;
        if q >= T::ONE {
            return T::ZERO;
        }
        let base = T::ONE - q;
        T::from_f64(30.0) / (T::PI * self.support * self.support * self.support) * base * base
    }

    /// Lage des Teilchens mit der Nummer `index`: zuerst die Flüssigkeit,
    /// dann die Wände.
    #[inline]
    fn point(&self, index: usize) -> Vec2<T> {
        let count = self.positions.len();
        if index < count { self.positions[index] } else { self.walls[index - count] }
    }

    /// Hält ein Teilchen im Becken und außerhalb der Hindernisse.
    fn confine(&self, mut point: Vec2<T>) -> Vec2<T> {
        let margin = self.spacing * T::from_f64(0.25);
        for axis in 0..2 {
            let (low, high) = (self.bounds.min[axis] + margin, self.bounds.max[axis] - margin);
            point[axis] = if point[axis].is_nan() { low } else { point[axis].max(low).min(high.max(low)) };
        }
        for block in &self.obstacles {
            let grown = block.expanded(margin);
            if !grown.contains(point) {
                continue;
            }
            // Auf die nächstgelegene Seite hinausschieben.
            let exits = [
                (point[0] - grown.min[0], 0, grown.min[0]),
                (grown.max[0] - point[0], 0, grown.max[0]),
                (point[1] - grown.min[1], 1, grown.min[1]),
                (grown.max[1] - point[1], 1, grown.max[1]),
            ];
            let mut best = exits[0];
            for exit in exits {
                if exit.0 < best.0 {
                    best = exit;
                }
            }
            point[best.1] = best.2;
        }
        point
    }

    /// Sortiert alle Teilchen in Zellen und sammelt zu jedem – auch zu den
    /// Wandteilchen – die Nachbarn in Reichweite.
    fn find_neighbours(&mut self) {
        let count = self.positions.len();
        let total = count + self.walls.len();
        let border = Vector::splat(self.spacing * T::from_f64(WALL_LAYERS));
        let origin = self.bounds.min - border;
        let extent = self.bounds.max - origin + border;
        let cells_along = |length: T| {
            let cells = (length / self.support).ceil().to_f64();
            if cells.is_finite() { (cells as usize).clamp(1, 4096) } else { 1 }
        };
        let (width, height) = (cells_along(extent[0]), cells_along(extent[1]));
        let cell_of = |point: Vec2<T>| {
            let along = |value: T, cells: usize| {
                let scaled = (value / self.support).floor().to_f64();
                if scaled.is_finite() { (scaled.max(0.0) as usize).min(cells - 1) } else { 0 }
            };
            (along(point[0] - origin[0], width), along(point[1] - origin[1], height))
        };

        // Zählsortierung nach Zellen.
        self.cell_start.clear();
        self.cell_start.resize(width * height + 1, 0);
        for index in 0..total {
            let (cx, cy) = cell_of(self.point(index));
            self.cell_start[cy * width + cx + 1] += 1;
        }
        for cell in 0..width * height {
            self.cell_start[cell + 1] += self.cell_start[cell];
        }
        let mut cursor = self.cell_start.clone();
        self.cell_items.clear();
        self.cell_items.resize(total, 0);
        for index in 0..total {
            let (cx, cy) = cell_of(self.point(index));
            let slot = &mut cursor[cy * width + cx];
            self.cell_items[*slot] = index;
            *slot += 1;
        }

        self.neighbour_start.clear();
        self.neighbours.clear();
        for index in 0..total {
            self.neighbour_start.push(self.neighbours.len());
            let here = self.point(index);
            let (cx, cy) = cell_of(here);
            for ny in cy.saturating_sub(1)..=(cy + 1).min(height - 1) {
                for nx in cx.saturating_sub(1)..=(cx + 1).min(width - 1) {
                    let cell = ny * width + nx;
                    for &other in &self.cell_items[self.cell_start[cell]..self.cell_start[cell + 1]] {
                        if other != index && (self.point(other) - here).norm() < self.support {
                            self.neighbours.push(other);
                        }
                    }
                }
            }
        }
        self.neighbour_start.push(self.neighbours.len());
    }

    /// Ein Durchlauf der Dichtekorrektur.
    fn correct_density(&mut self) {
        let count = self.positions.len();
        let softness = self.params.relaxation.max(T::ZERO) / (self.support * self.support);
        // Dichte und Druckfaktor an jedem Teilchen, auch an denen der Wand.
        // Eine Wand, die nur den Druck ihres Gegenübers spiegelt, trägt kein
        // Druckgefälle mit der Tiefe: Der Spalte an der Wand fehlt dann die
        // halbe Stütze gegen die Schwerkraft, sie sinkt, und es entsteht
        // eine Umwälzung, die nie zur Ruhe kommt. Mit eigener Dichte wächst
        // der Druck in der Wand mit der Tiefe wie in der Flüssigkeit.
        for index in 0..count + self.walls.len() {
            let here = self.point(index);
            let movable = index < count;
            let mut density = self.density_kernel(T::ZERO);
            let mut own_gradient: Vec2<T> = Vector::ZERO;
            let mut squares = T::ZERO;
            let mut wet = movable;
            for &other in &self.neighbours[self.neighbour_start[index]..self.neighbour_start[index + 1]] {
                let offset = here - self.point(other);
                let distance = offset.norm();
                density += self.density_kernel(distance);
                wet |= other < count;
                if distance > T::ZERO {
                    let slope = self.gradient_kernel(distance) / self.rest_density;
                    own_gradient -= offset * (slope / distance);
                    squares += slope * slope;
                }
            }
            let relative = density / self.rest_density;
            if movable {
                self.densities[index] = relative;
            }
            // Nur Druck, kein Zug: Zu dünn besetzte Stellen (die freie
            // Oberfläche) ziehen nichts an. Wandteilchen ohne Flüssigkeit in
            // Reichweite drücken auf nichts.
            let excess = if wet { (relative - T::ONE).max(T::ZERO) } else { T::ZERO };
            self.lambdas[index] = -excess / (squares + own_gradient.norm_sqr() + softness);
        }

        let reference = self.density_kernel(self.support * T::from_f64(0.3));
        let clump = self.params.anti_clumping.max(T::ZERO) * self.support * self.support;
        for index in 0..count {
            let here = self.positions[index];
            let mut delta: Vec2<T> = Vector::ZERO;
            for &other in &self.neighbours[self.neighbour_start[index]..self.neighbour_start[index + 1]] {
                let offset = here - self.point(other);
                let distance = offset.norm();
                // Deckungsgleiche Teilchen bekommen eine feste Richtung,
                // sonst kämen sie nie wieder auseinander.
                let direction = if distance > T::ZERO {
                    offset * distance.recip()
                } else if index < other {
                    Vector([T::ONE, T::ZERO])
                } else {
                    Vector([-T::ONE, T::ZERO])
                };
                let partner = self.lambdas[other];
                let ratio = self.density_kernel(distance) / reference;
                let ratio = ratio * ratio;
                let weight = self.lambdas[index] + partner - clump * ratio * ratio;
                delta -= direction * (weight * self.gradient_kernel(distance) / self.rest_density);
            }
            self.deltas[index] = delta;
        }
        for index in 0..count {
            // Höchstens ein Teilchenabstand je Durchlauf: hält die Korrektur
            // auch bei extremer Verdichtung im Rahmen.
            let mut delta = self.deltas[index];
            let length = delta.norm();
            if length > self.spacing {
                delta *= self.spacing / length;
            } else if length.is_nan() {
                delta = Vector::ZERO;
            }
            self.positions[index] = self.confine(self.positions[index] + delta);
        }
    }

    /// Gleicht die Geschwindigkeit an die der Nachbarn an (XSPH).
    fn smooth_velocities(&mut self) {
        let count = self.positions.len();
        let strength = self.params.viscosity.max(T::ZERO).min(T::ONE);
        if strength == T::ZERO {
            return;
        }
        for index in 0..count {
            let here = self.positions[index];
            let mut change: Vec2<T> = Vector::ZERO;
            for &other in &self.neighbours[self.neighbour_start[index]..self.neighbour_start[index + 1]] {
                if other < count {
                    let weight = self.density_kernel((here - self.positions[other]).norm()) / self.rest_density;
                    change += (self.velocities[other] - self.velocities[index]) * weight;
                }
            }
            self.deltas[index] = change * strength;
        }
        for index in 0..count {
            self.velocities[index] += self.deltas[index];
        }
    }

    /// Rechnet um `dt` weiter. Ein nicht positiver oder nicht endlicher
    /// Zeitschritt ändert nichts.
    pub fn step(&mut self, dt: T) {
        if dt <= T::ZERO || !dt.is_finite() || self.positions.is_empty() {
            return;
        }
        let substeps = self.params.substeps.max(1);
        let h = dt / T::from_f64(f64::from(substeps));
        for _ in 0..substeps {
            for index in 0..self.positions.len() {
                self.velocities[index] += self.params.gravity * h;
                self.previous[index] = self.positions[index];
                self.positions[index] = self.confine(self.positions[index] + self.velocities[index] * h);
            }
            self.lambdas.resize(self.positions.len() + self.walls.len(), T::ZERO);
            self.find_neighbours();
            for _ in 0..self.params.iterations.max(1) {
                self.correct_density();
            }
            for index in 0..self.positions.len() {
                self.velocities[index] = (self.positions[index] - self.previous[index]) * h.recip();
            }
            self.smooth_velocities();
        }
    }
}


#[cfg(all(test, any(feature = "std", feature = "libm")))]
mod confinement_tests {
    use super::*;

    /// Ein Hindernis steht auf dem Beckenboden. Die kürzeste geometrische
    /// Flucht nach unten ist unzulässig: Das Ergebnis muss gleichzeitig
    /// innerhalb des Beckens und außerhalb des Hindernisses liegen.
    ///
    /// Direkter Test von confine, unabhängig von Dichte, Druck, Zeitschritt
    /// und Darstellung. Der Eingang darf im Hindernis liegen, denn genau
    /// solche Lagen muss die Kollisionskorrektur auflösen.
    fn floor_obstacle_keeps_particles_in_tank<T: RealElementary>() {
        let f = T::from_f64;
        let mut fluid = ParticleFluid::<T>::new(
            Aabb { min: Vector([f(0.0), f(0.0)]), max: Vector([f(4.0), f(2.4)]) },
            f(0.04),
        );
        fluid.add_obstacle(Vector([f(2.2), f(0.0)]), Vector([f(2.4), f(0.4)]));
        let margin = fluid.spacing() * f(0.25);
        let low = fluid.bounds().min + Vector::splat(margin);
        let high = fluid.bounds().max - Vector::splat(margin);
        let tolerance = T::EPSILON * f(16.0);

        // Zuerst das konkrete Gegenbeispiel, dann Punkte nahe beiden Ecken.
        for x in [2.3, 2.21, 2.39] {
            let input = Vector([f(x), margin]);
            let output = fluid.confine(input);
            assert!(output.0.iter().all(|v| v.is_finite()));
            assert!(fluid.bounds().contains(output),
                "confine schiebt Teilchen unter den Boden: Eingang {:?}, Ergebnis {:?}",
                input.0.map(|v| v.to_f64()), output.0.map(|v| v.to_f64()));
            for axis in 0..2 {
                assert!(output[axis] >= low[axis] - tolerance
                    && output[axis] <= high[axis] + tolerance,
                    "Mindestabstand zum Beckenrand verletzt, Achse {axis}");
            }
            assert!(!fluid.obstacles()[0].contains(output),
                "Teilchen blieb im Hindernis");
            let again = fluid.confine(output);
            assert!((again - output).norm() <= tolerance,
                "Wiederholtes confine verschiebt einen bereits gültigen Kontakt");
        }
    }

    #[test]
    fn floor_obstacle_keeps_particles_in_tank_f32() {
        floor_obstacle_keeps_particles_in_tank::<f32>();
    }

    #[test]
    fn floor_obstacle_keeps_particles_in_tank_f64() {
        floor_obstacle_keeps_particles_in_tank::<f64>();
    }
}
