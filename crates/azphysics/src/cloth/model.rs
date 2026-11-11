//! Das Tuch: Teilchen, Verbindungen, Dreiecke und Einstellungen.

use crate::geometry::Vec3;
use alloc::vec;
use alloc::vec::Vec;
use azmath::Vector;
use azmath::scalar::RealElementary;

/// Einstellungen der Simulation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClothParams<T> {
    /// Schwerebeschleunigung.
    pub gravity: Vec3<T>,
    /// Teilschritte je Zeitschritt. Mehr Teilschritte machen das Tuch
    /// steifer und die Kollisionen genauer.
    pub substeps: u32,
    /// Durchläufe über alle Abstandsbedingungen je Teilschritt.
    pub iterations: u32,
    /// Nachgiebigkeit gegen Dehnen; null ist undehnbar.
    pub stretch_compliance: T,
    /// Nachgiebigkeit gegen Biegen; größere Werte ergeben weicheren Fall.
    pub bend_compliance: T,
    /// Dämpfung der Geschwindigkeit je Sekunde (0 = keine).
    pub damping: T,
    /// Dicke des Tuchs: Mindestabstand zu Hindernissen und zu sich selbst.
    pub thickness: T,
    /// Reibung an Hindernissen, 0 (glatt) bis 1 (haftend).
    pub friction: T,
    /// Abstoßung zwischen Tuch und Dreiecken (auch den eigenen).
    pub self_collision: bool,
    /// Sicherung gegen Durchdringung mit kontinuierlicher
    /// Kollisionserkennung.
    pub ccd: bool,
}

impl<T: RealElementary> Default for ClothParams<T> {
    fn default() -> Self {
        Self {
            gravity: Vector([T::ZERO, T::from_f64(-9.81), T::ZERO]),
            substeps: 8,
            iterations: 2,
            stretch_compliance: T::ZERO,
            bend_compliance: T::from_f64(1.0e-3),
            damping: T::from_f64(0.3),
            thickness: T::from_f64(0.01),
            friction: T::from_f64(0.3),
            self_collision: true,
            ccd: true,
        }
    }
}

/// Einfaches Hindernis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Collider<T> {
    /// Halbraum: alles auf der Seite, in die `normal` zeigt, ist frei.
    Plane {
        /// Ein Punkt der Ebene.
        point: Vec3<T>,
        /// Normale der Länge eins, zeigt in den freien Raum.
        normal: Vec3<T>,
    },
    /// Kugel.
    Sphere {
        /// Mittelpunkt.
        center: Vec3<T>,
        /// Radius.
        radius: T,
    },
}

/// Was in einem Zeitschritt geschah.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StepStats {
    /// Abstoßungen zwischen Teilchen und Dreiecken.
    pub repulsions: u32,
    /// Von der kontinuierlichen Prüfung gefundene Berührungen.
    pub ccd_hits: u32,
    /// Teilschritte, in denen alle Teilchen zurückgesetzt werden mussten.
    pub full_rollbacks: u32,
}

/// Abstandsbedingung zwischen zwei Teilchen.
#[derive(Clone, Copy, Debug)]
pub(super) struct Link<T> {
    pub(super) i: usize,
    pub(super) j: usize,
    pub(super) rest: T,
    /// Biegebedingung (überspringt ein Teilchen) statt Dehnbedingung.
    pub(super) bend: bool,
}

/// Ein Tuch aus Teilchen.
#[derive(Clone, Debug)]
pub struct Cloth<T> {
    pub(super) positions: Vec<Vec3<T>>,
    pub(super) previous: Vec<Vec3<T>>,
    pub(super) velocities: Vec<Vec3<T>>,
    pub(super) inverse_mass: Vec<T>,
    pub(super) links: Vec<Link<T>>,
    pub(super) lambdas: Vec<T>,
    pub(super) triangles: Vec<[usize; 3]>,
    pub(super) edges: Vec<[usize; 2]>,
    cloth_triangles: usize,
    columns: usize,
    /// Einstellungen; dürfen zwischen den Schritten geändert werden.
    pub params: ClothParams<T>,
}

impl<T: RealElementary> Cloth<T> {
    /// Rechteckiges Tuch aus `columns × rows` Teilchen (jeweils mindestens
    /// zwei).
    ///
    /// Teilchen `(spalte, zeile)` liegt bei
    /// `origin + spalte·across + zeile·down`; `across` und `down` sind also
    /// die Abstände zum Nachbarn. Jedes Teilchen hat die Masse
    /// `particle_mass`.
    pub fn grid(
        columns: usize,
        rows: usize,
        origin: Vec3<T>,
        across: Vec3<T>,
        down: Vec3<T>,
        particle_mass: T,
    ) -> Self {
        let (columns, rows) = (columns.max(2), rows.max(2));
        let count = columns * rows;
        let positions: Vec<Vec3<T>> = (0..count)
            .map(|index| {
                let (column, row) = (index % columns, index / columns);
                origin + across * T::from_f64(column as f64) + down * T::from_f64(row as f64)
            })
            .collect();
        let at = |column: usize, row: usize| row * columns + column;

        let mut links = Vec::new();
        let mut link = |i: usize, j: usize, bend: bool| {
            links.push(Link { i, j, rest: (positions[i] - positions[j]).norm(), bend });
        };
        let mut triangles = Vec::new();
        for row in 0..rows {
            for column in 0..columns {
                // Dehnung: Nachbarn längs und quer.
                if column + 1 < columns {
                    link(at(column, row), at(column + 1, row), false);
                }
                if row + 1 < rows {
                    link(at(column, row), at(column, row + 1), false);
                }
                // Scherung: beide Diagonalen der Zelle.
                if column + 1 < columns && row + 1 < rows {
                    link(at(column, row), at(column + 1, row + 1), false);
                    link(at(column + 1, row), at(column, row + 1), false);
                    triangles.push([at(column, row), at(column + 1, row), at(column, row + 1)]);
                    triangles.push([at(column + 1, row), at(column + 1, row + 1), at(column, row + 1)]);
                }
                // Biegung: übernächste Nachbarn.
                if column + 2 < columns {
                    link(at(column, row), at(column + 2, row), true);
                }
                if row + 2 < rows {
                    link(at(column, row), at(column, row + 2), true);
                }
            }
        }

        let inverse = if particle_mass > T::ZERO { particle_mass.recip() } else { T::ZERO };
        let mut cloth = Self {
            previous: positions.clone(),
            velocities: vec![Vector::ZERO; count],
            inverse_mass: vec![inverse; count],
            lambdas: vec![T::ZERO; links.len()],
            cloth_triangles: triangles.len(),
            positions,
            links,
            triangles,
            edges: Vec::new(),
            columns,
            params: ClothParams::default(),
        };
        cloth.rebuild_edges();
        cloth
    }

    /// Sammelt die Kanten aller Dreiecke, jede einmal.
    fn rebuild_edges(&mut self) {
        self.edges.clear();
        for triangle in &self.triangles {
            for (a, b) in [(0, 1), (1, 2), (2, 0)] {
                let (i, j) = (triangle[a], triangle[b]);
                self.edges.push([i.min(j), i.max(j)]);
            }
        }
        self.edges.sort_unstable();
        self.edges.dedup();
    }

    /// Nummer des Teilchens in Spalte `column` und Zeile `row` eines mit
    /// [`grid`](Self::grid) erzeugten Tuchs.
    #[inline]
    pub fn index(&self, column: usize, row: usize) -> usize {
        row * self.columns + column
    }

    /// Lagen aller Teilchen (Tuch zuerst, dann die Ecken der Hindernisse).
    #[inline]
    pub fn positions(&self) -> &[Vec3<T>] {
        &self.positions
    }

    /// Geschwindigkeiten aller Teilchen.
    #[inline]
    pub fn velocities(&self) -> &[Vec3<T>] {
        &self.velocities
    }

    /// Alle Dreiecke als Teilchennummern: zuerst die des Tuchs, dann die der
    /// Hindernisse.
    #[inline]
    pub fn triangles(&self) -> &[[usize; 3]] {
        &self.triangles
    }

    /// Anzahl der Dreiecke, die zum Tuch gehören.
    #[inline]
    pub fn cloth_triangle_count(&self) -> usize {
        self.cloth_triangles
    }

    /// Hält ein Teilchen fest: Es bewegt sich nicht mehr von selbst.
    pub fn pin(&mut self, index: usize) {
        if let Some(inverse) = self.inverse_mass.get_mut(index) {
            *inverse = T::ZERO;
        }
        if let Some(velocity) = self.velocities.get_mut(index) {
            *velocity = Vector::ZERO;
        }
    }

    /// Versetzt ein Teilchen, etwa um ein festgehaltenes zu führen.
    pub fn set_position(&mut self, index: usize, position: Vec3<T>) {
        if let Some(slot) = self.positions.get_mut(index) {
            *slot = position;
        }
    }

    /// Setzt die Geschwindigkeit eines Teilchens.
    pub fn set_velocity(&mut self, index: usize, velocity: Vec3<T>) {
        if let (Some(slot), Some(&inverse)) = (self.velocities.get_mut(index), self.inverse_mass.get(index)) {
            if inverse > T::ZERO {
                *slot = velocity;
            }
        }
    }

    /// Gibt allen beweglichen Teilchen dieselbe Geschwindigkeit.
    pub fn set_all_velocities(&mut self, velocity: Vec3<T>) {
        for index in 0..self.velocities.len() {
            self.set_velocity(index, velocity);
        }
    }

    /// Fügt ein ruhendes Dreieck als Hindernis hinzu.
    ///
    /// Das Tuch hält Abstand zu ihm und kann es – mit eingeschalteter
    /// Sicherung – nicht durchdringen.
    pub fn add_obstacle_triangle(&mut self, a: Vec3<T>, b: Vec3<T>, c: Vec3<T>) {
        let first = self.positions.len();
        for corner in [a, b, c] {
            self.positions.push(corner);
            self.previous.push(corner);
            self.velocities.push(Vector::ZERO);
            self.inverse_mass.push(T::ZERO);
        }
        self.triangles.push([first, first + 1, first + 2]);
        self.rebuild_edges();
    }

    /// Größte relative Längenänderung einer Dehnverbindung gegenüber der
    /// Ruhelänge; null heißt ungedehnt.
    pub fn max_strain(&self) -> T {
        self.links.iter().filter(|link| !link.bend && link.rest > T::ZERO).fold(
            T::ZERO,
            |worst, link| {
                let length = (self.positions[link.i] - self.positions[link.j]).norm();
                worst.max(((length - link.rest) / link.rest).abs())
            },
        )
    }
}
