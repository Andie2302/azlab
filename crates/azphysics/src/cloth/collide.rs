//! Kollision des Tuchs mit Dreiecken: Abstoßung und Sicherung gegen
//! Durchdringung.

use super::model::{Cloth, StepStats};
use crate::broadphase::SpatialHash;
use crate::ccd::{edge_edge, vertex_triangle};
use crate::geometry::{Aabb, Vec3, closest_point_on_triangle, segment_intersects_triangle, triangle_normal};
use alloc::vec::Vec;
use azmath::scalar::RealElementary;

/// Durchläufe der Sicherung, bevor alles zurückgesetzt wird.
const MAX_PASSES: usize = 8;

impl<T: RealElementary> Cloth<T> {
    /// Kantenlänge der Suchzellen: etwa eine Maschenweite.
    fn cell_size(&self) -> T {
        let (sum, count) = self.links.iter().filter(|link| !link.bend).fold(
            (T::ZERO, T::ZERO),
            |(sum, count), link| (sum + link.rest, count + T::ONE),
        );
        let mesh = if count > T::ZERO { sum / count } else { T::ONE };
        mesh.max(self.params.thickness * T::from_f64(2.0))
    }

    /// Hüllquader von Punkten in Anfangs- und Endlage des Teilschritts.
    fn swept(&self, indices: &[usize], margin: T) -> Aabb<T, 3> {
        let points = indices.iter().flat_map(|&i| [self.previous[i], self.positions[i]]);
        Aabb::from_points(points)
            .unwrap_or(Aabb { min: self.positions[0], max: self.positions[0] })
            .expanded(margin)
    }

    /// `true`, wenn sich keins der Teilchen von selbst bewegen kann.
    fn all_fixed(&self, indices: &[usize]) -> bool {
        indices.iter().all(|&i| self.inverse_mass[i] <= T::ZERO)
    }

    /// Schiebt Teilchen und Dreiecke auseinander, die sich näher kommen als
    /// die Dicke des Tuchs.
    pub(super) fn repel(&mut self, stats: &mut StepStats) {
        let thickness = self.params.thickness;
        let hash = SpatialHash::build(
            self.cell_size(),
            self.triangles.iter().enumerate().map(|(id, tri)| (id, self.swept(tri, thickness))),
        );
        let mut candidates = Vec::new();
        for vertex in 0..self.positions.len() {
            hash.query(self.swept(&[vertex], thickness), &mut candidates);
            for &id in &candidates {
                let triangle = self.triangles[id];
                if triangle.contains(&vertex) || (self.all_fixed(&[vertex]) && self.all_fixed(&triangle)) {
                    continue;
                }
                let [a, b, c] = triangle;
                let point = self.positions[vertex];
                let (closest, weights) =
                    closest_point_on_triangle(point, self.positions[a], self.positions[b], self.positions[c]);
                let offset = point - closest;
                let distance = offset.norm();
                if distance >= thickness {
                    continue;
                }
                // Richtung: weg vom Dreieck, auf die Seite, auf der das
                // Teilchen am Anfang des Teilschritts lag.
                let before = self.previous[vertex]
                    - closest_point_on_triangle(
                        self.previous[vertex],
                        self.previous[a],
                        self.previous[b],
                        self.previous[c],
                    )
                    .0;
                let normal = triangle_normal(self.positions[a], self.positions[b], self.positions[c]);
                let side = if normal.dot(before) < T::ZERO { -normal } else { normal };
                let (direction, gap) = match (offset.normalize(), side.normalize()) {
                    (Some(away), _) if away.dot(side) >= T::ZERO => (away, thickness - distance),
                    (_, Some(side)) => (side, thickness + distance),
                    _ => continue,
                };
                let inverse = [
                    self.inverse_mass[vertex],
                    self.inverse_mass[a],
                    self.inverse_mass[b],
                    self.inverse_mass[c],
                ];
                let share = inverse[0]
                    + weights[0] * weights[0] * inverse[1]
                    + weights[1] * weights[1] * inverse[2]
                    + weights[2] * weights[2] * inverse[3];
                if share <= T::ZERO {
                    continue;
                }
                let push = direction * (gap / share);
                self.positions[vertex] += push * inverse[0];
                self.positions[a] -= push * (weights[0] * inverse[1]);
                self.positions[b] -= push * (weights[1] * inverse[2]);
                self.positions[c] -= push * (weights[2] * inverse[3]);
                stats.repulsions += 1;
            }
        }
    }

    /// Setzt Teilchen auf ihre Lage vom Anfang des Teilschritts zurück.
    fn roll_back(&mut self, indices: &[usize]) {
        for &index in indices {
            self.positions[index] = self.previous[index];
        }
    }

    /// Bewegung eines Teilchens im Teilschritt.
    #[inline]
    fn motion(&self, index: usize) -> [Vec3<T>; 2] {
        [self.previous[index], self.positions[index]]
    }

    /// Sicherung gegen Durchdringung.
    ///
    /// Prüft jede Bewegung Punkt gegen Dreieck und Kante gegen Kante. Wo
    /// sich Teile im Teilschritt berühren würden, werden die beteiligten
    /// Teilchen zurückgesetzt. Das wird wiederholt, weil ein Zurücksetzen
    /// die Bewegungen der Nachbarn ändert; bleibt es dabei nicht ruhig,
    /// werden alle Teilchen zurückgesetzt. Ohne Bewegung gibt es keine neue
    /// Durchdringung.
    pub(super) fn prevent_penetration(&mut self, stats: &mut StepStats) {
        let tolerance = self.params.thickness * T::from_f64(0.1);
        let cell = self.cell_size();
        let triangle_hash = SpatialHash::build(
            cell,
            self.triangles.iter().enumerate().map(|(id, tri)| (id, self.swept(tri, tolerance))),
        );
        let edge_hash = SpatialHash::build(
            cell,
            self.edges.iter().enumerate().map(|(id, edge)| (id, self.swept(edge, tolerance))),
        );
        let mut candidates = Vec::new();

        for _ in 0..MAX_PASSES {
            let mut clean = true;
            for vertex in 0..self.positions.len() {
                triangle_hash.query(self.swept(&[vertex], tolerance), &mut candidates);
                for &id in &candidates {
                    let triangle = self.triangles[id];
                    if triangle.contains(&vertex)
                        || (self.all_fixed(&[vertex]) && self.all_fixed(&triangle))
                    {
                        continue;
                    }
                    let [a, b, c] = triangle;
                    let involved = [vertex, a, b, c];
                    if involved.iter().all(|&i| self.positions[i] == self.previous[i]) {
                        continue;
                    }
                    let hit = vertex_triangle(
                        self.motion(vertex),
                        self.motion(a),
                        self.motion(b),
                        self.motion(c),
                        tolerance,
                    );
                    if hit.is_some() {
                        self.roll_back(&involved);
                        stats.ccd_hits += 1;
                        clean = false;
                    }
                }
            }
            for first in 0..self.edges.len() {
                let [a, b] = self.edges[first];
                edge_hash.query(self.swept(&[a, b], tolerance), &mut candidates);
                for &second in &candidates {
                    let [c, d] = self.edges[second];
                    // Jedes Paar einmal; Kanten mit gemeinsamem Teilchen
                    // können sich nicht durchdringen.
                    if second <= first || a == c || a == d || b == c || b == d {
                        continue;
                    }
                    let involved = [a, b, c, d];
                    if self.all_fixed(&involved)
                        || involved.iter().all(|&i| self.positions[i] == self.previous[i])
                    {
                        continue;
                    }
                    let hit = edge_edge(
                        self.motion(a),
                        self.motion(b),
                        self.motion(c),
                        self.motion(d),
                        tolerance,
                    );
                    if hit.is_some() {
                        self.roll_back(&involved);
                        stats.ccd_hits += 1;
                        clean = false;
                    }
                }
            }
            if clean {
                return;
            }
        }
        // Nicht zur Ruhe gekommen: alles zurücksetzen.
        let all: Vec<usize> = (0..self.positions.len()).collect();
        self.roll_back(&all);
        stats.full_rollbacks += 1;
    }

    /// Zählt die Stellen, an denen eine Kante ein fremdes Dreieck
    /// durchstößt.
    ///
    /// Null heißt: Das Netz durchdringt weder sich selbst noch ein
    /// Hindernis. Gedacht zum Nachprüfen, nicht für jeden Zeitschritt.
    pub fn intersections(&self) -> usize {
        let hash = SpatialHash::build(
            self.cell_size(),
            self.triangles.iter().enumerate().map(|(id, tri)| {
                let points = tri.iter().map(|&i| self.positions[i]);
                let bounds = Aabb::from_points(points)
                    .unwrap_or(Aabb { min: self.positions[0], max: self.positions[0] });
                (id, bounds)
            }),
        );
        let mut candidates = Vec::new();
        let mut count = 0;
        for &[p, q] in &self.edges {
            let Some(bounds) = Aabb::from_points([self.positions[p], self.positions[q]]) else {
                continue;
            };
            hash.query(bounds, &mut candidates);
            for &id in &candidates {
                let [a, b, c] = self.triangles[id];
                if [a, b, c].contains(&p) || [a, b, c].contains(&q) {
                    continue;
                }
                if segment_intersects_triangle(
                    self.positions[p],
                    self.positions[q],
                    self.positions[a],
                    self.positions[b],
                    self.positions[c],
                ) {
                    count += 1;
                }
            }
        }
        count
    }
}
