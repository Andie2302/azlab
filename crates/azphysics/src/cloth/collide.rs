//! Kollision des Tuchs mit Dreiecken: Abstoßung und Sicherung gegen
//! Durchdringung.

use super::model::{Cloth, StepStats, StepPhase};
use crate::broadphase::SpatialHash;
use crate::ccd::{edge_edge, vertex_triangle};
use crate::geometry::{
    Aabb, Vec3, closest_point_on_triangle, point_triangle_distance, segment_intersects_triangle,
    segment_segment_distance, triangle_normal,
};
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

    /// Sorgt dafür, dass die Nachbarschaftslisten zum aktuellen Stand
    /// passen.
    ///
    /// Die Listen nennen alle Paare aus Teilchen und Dreieck sowie aus zwei
    /// Kanten, die sich nahe genug sind, um sich berühren zu können. Sie
    /// werden mit einem Sicherheitsabstand gebaut und bleiben gültig,
    /// solange sich kein Teilchen weiter als den halben Sicherheitsabstand
    /// aus seinem Bereich beim Bau entfernt hat; erst dann wird neu gebaut.
    /// Ein Paar, das nicht in der Liste steht, kann sich also nicht
    /// berühren.
    fn refresh_proximity(&mut self) {
        let thickness = self.params.thickness;
        let near = &self.proximity;
        let current = near.valid
            && near.thickness == thickness
            && near.anchors.len() == self.positions.len()
            && near
                .anchors
                .iter()
                .zip(&self.positions)
                .all(|(anchor, &position)| anchor.expanded(near.slack).contains(position));
        if current {
            return;
        }

        let cell = self.cell_size();
        let slack = cell * T::from_f64(0.25);
        let reach = slack + thickness;
        let mut near = core::mem::take(&mut self.proximity);
        near.anchors.clear();
        near.anchors.extend(
            self.previous.iter().zip(&self.positions).map(|(&from, &to)| Aabb {
                min: from.zip_map(to, |a, b| a.min(b)),
                max: from.zip_map(to, |a, b| a.max(b)),
            }),
        );
        let hull = |indices: &[usize]| {
            indices[1..].iter().fold(near.anchors[indices[0]], |hull, &i| hull.union(near.anchors[i])).expanded(reach)
        };

        let triangle_boxes: Vec<Aabb<T, 3>> = self.triangles.iter().map(|tri| hull(tri)).collect();
        let mut hash = SpatialHash::build(cell, triangle_boxes.iter().copied().enumerate());
        let mut candidates = Vec::new();
        near.vertex_triangle.clear();
        for vertex in 0..self.positions.len() {
            let bounds = hull(&[vertex]);
            hash.query(bounds, &mut candidates);
            for &id in &candidates {
                let triangle = self.triangles[id];
                if !triangle.contains(&vertex)
                    && triangle_boxes[id].overlaps(bounds)
                    && !(self.all_fixed(&[vertex]) && self.all_fixed(&triangle))
                {
                    near.vertex_triangle.push((vertex, id));
                }
            }
        }

        let edge_boxes: Vec<Aabb<T, 3>> = self.edges.iter().map(|edge| hull(edge)).collect();
        let mut hash = SpatialHash::build(cell, edge_boxes.iter().copied().enumerate());
        near.edge_edge.clear();
        for first in 0..self.edges.len() {
            let [a, b] = self.edges[first];
            hash.query(edge_boxes[first], &mut candidates);
            for &second in &candidates {
                let [c, d] = self.edges[second];
                // Jedes Paar einmal; Kanten mit gemeinsamem Teilchen können
                // sich nicht durchdringen.
                if second > first
                    && a != c
                    && a != d
                    && b != c
                    && b != d
                    && edge_boxes[first].overlaps(edge_boxes[second])
                    && !self.all_fixed(&[a, b, c, d])
                {
                    near.edge_edge.push((first, second));
                }
            }
        }
        near.slack = slack;
        near.thickness = thickness;
        near.valid = true;
        self.proximity = near;
    }

    /// Schiebt Teilchen und Dreiecke auseinander, die sich näher kommen als
    /// die Dicke des Tuchs.
    pub(super) fn repel(&mut self, stats: &mut StepStats, observe: &mut impl FnMut(StepPhase)) {
        let thickness = self.params.thickness;
        self.refresh_proximity();
        observe(StepPhase::Proximity);
        let near = core::mem::take(&mut self.proximity);
        for &(vertex, id) in &near.vertex_triangle {
            let [a, b, c] = self.triangles[id];
            let point = self.positions[vertex];
            // Schnelle Vorprüfung: außerhalb des Hüllquaders ist das
            // Teilchen weiter weg als die Dicke.
            let bounds = Aabb {
                min: self.positions[a].zip_map(self.positions[b], |x, y| x.min(y)),
                max: self.positions[a].zip_map(self.positions[b], |x, y| x.max(y)),
            };
            let bounds = Aabb {
                min: bounds.min.zip_map(self.positions[c], |x, y| x.min(y)),
                max: bounds.max.zip_map(self.positions[c], |x, y| x.max(y)),
            };
            if !bounds.expanded(thickness).contains(point) {
                continue;
            }
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
                - closest_point_on_triangle(self.previous[vertex], self.previous[a], self.previous[b], self.previous[c])
                    .0;
            let normal = triangle_normal(self.positions[a], self.positions[b], self.positions[c]);
            let side = if normal.dot(before) < T::ZERO { -normal } else { normal };
            let (direction, gap) = match (offset.normalize(), side.normalize()) {
                (Some(away), _) if away.dot(side) >= T::ZERO => (away, thickness - distance),
                (_, Some(side)) => (side, thickness + distance),
                _ => continue,
            };
            let inverse =
                [self.inverse_mass[vertex], self.inverse_mass[a], self.inverse_mass[b], self.inverse_mass[c]];
            let share = inverse[0]
                + weights[0] * weights[0] * inverse[1]
                + weights[1] * weights[1] * inverse[2]
                + weights[2] * weights[2] * inverse[3];
            if share <= T::ZERO {
                continue;
            }
            let push = direction * (gap / share);
            let moves = [
                (vertex, push * inverse[0]),
                (a, -push * (weights[0] * inverse[1])),
                (b, -push * (weights[1] * inverse[2])),
                (c, -push * (weights[2] * inverse[3])),
            ];
            for (index, shift) in moves {
                self.positions[index] += shift;
                self.pushes[index] += shift;
            }
            stats.repulsions += 1;
        }
        self.proximity = near;
    }

    /// Setzt Teilchen auf ihre Lage vom Anfang des Teilschritts zurück.
    fn roll_back(&mut self, indices: &[usize]) {
        for &index in indices {
            self.positions[index] = self.previous[index];
        }
    }

    /// Kürzt die Bewegung der beteiligten Teilchen so, dass sie kurz vor der
    /// Berührung zur Zeit `contact` enden.
    ///
    /// Die Teilchen bleiben auf ihrer Bahn und halten dort an, wo `gap`
    /// (der Abstand der beiden Teile in der neuen Lage) noch deutlich über
    /// der Toleranz liegt. Gelingt das nicht, oder mit `gentle = false`,
    /// kehren sie an den Anfang des Teilschritts zurück.
    fn retreat(&mut self, involved: &[usize; 4], contact: T, tolerance: T, gentle: bool, gap: impl Fn(&Self) -> T) {
        if gentle {
            let ends = involved.map(|index| self.positions[index]);
            for share in [0.8, 0.4] {
                let time = contact * T::from_f64(share);
                if time <= T::ZERO || time.is_nan() {
                    break;
                }
                for (&index, &end) in involved.iter().zip(&ends) {
                    self.positions[index] = self.previous[index] + (end - self.previous[index]) * time;
                }
                if gap(self) > tolerance * T::from_f64(3.0) {
                    return;
                }
            }
        }
        self.roll_back(involved);
    }

    /// Bewegung eines Teilchens im Teilschritt.
    #[inline]
    fn motion(&self, index: usize) -> [Vec3<T>; 2] {
        [self.previous[index], self.positions[index]]
    }

    /// Sicherung gegen Durchdringung.
    ///
    /// Prüft jede Bewegung Punkt gegen Dreieck und Kante gegen Kante. Wo
    /// sich Teile im Teilschritt berühren würden, wird die Bewegung der
    /// beteiligten Teilchen gekürzt: Sie halten kurz vor der Berührung an
    /// oder kehren an den Anfang des Teilschritts zurück. Das wird
    /// wiederholt, weil jede Kürzung die Bewegungen der Nachbarn ändert;
    /// bleibt es dabei nicht ruhig, werden alle Teilchen zurückgesetzt.
    /// Bewegungen werden nur gekürzt, nie verlängert, und der letzte
    /// Durchlauf hat alle Paare frei gefunden: Es gibt keine neue
    /// Durchdringung.
    pub(super) fn prevent_penetration(&mut self, stats: &mut StepStats, observe: &mut impl FnMut(StepPhase)) {
        let tolerance = self.params.thickness * T::from_f64(0.1);
        self.refresh_proximity();
        observe(StepPhase::Proximity);
        let near = core::mem::take(&mut self.proximity);
        // Hüllquader der Bewegungen, einmal berechnet. Zurücksetzen
        // verkleinert Bewegungen nur; die Quader bleiben also gültig. Paare
        // mit getrennten Quadern können sich nicht berühren.
        let paths: Vec<Aabb<T, 3>> = (0..self.positions.len()).map(|i| self.swept(&[i], tolerance)).collect();
        let hull = |indices: &[usize]| indices[1..].iter().fold(paths[indices[0]], |hull, &i| hull.union(paths[i]));
        let triangle_boxes: Vec<Aabb<T, 3>> = self.triangles.iter().map(|tri| hull(tri)).collect();
        let edge_boxes: Vec<Aabb<T, 3>> = self.edges.iter().map(|edge| hull(edge)).collect();

        let mut settled = false;
        for pass in 0..MAX_PASSES {
            // Erst behutsam kürzen, in späteren Durchläufen hart
            // zurücksetzen: Das kommt sicher zur Ruhe.
            let gentle = pass < MAX_PASSES / 2;
            let mut clean = true;
            for &(vertex, id) in &near.vertex_triangle {
                if !triangle_boxes[id].overlaps(paths[vertex]) {
                    continue;
                }
                let [a, b, c] = self.triangles[id];
                let involved = [vertex, a, b, c];
                if involved.iter().all(|&i| self.positions[i] == self.previous[i]) {
                    continue;
                }
                let hit =
                    vertex_triangle(self.motion(vertex), self.motion(a), self.motion(b), self.motion(c), tolerance);
                if let Some(contact) = hit {
                    self.retreat(&involved, contact, tolerance, gentle, |cloth| {
                        let at = |index: usize| cloth.positions[index];
                        point_triangle_distance(at(vertex), at(a), at(b), at(c))
                    });
                    stats.ccd_hits += 1;
                    clean = false;
                }
            }
            for &(first, second) in &near.edge_edge {
                if !edge_boxes[first].overlaps(edge_boxes[second]) {
                    continue;
                }
                let ([a, b], [c, d]) = (self.edges[first], self.edges[second]);
                let involved = [a, b, c, d];
                if involved.iter().all(|&i| self.positions[i] == self.previous[i]) {
                    continue;
                }
                let hit = edge_edge(self.motion(a), self.motion(b), self.motion(c), self.motion(d), tolerance);
                if let Some(contact) = hit {
                    self.retreat(&involved, contact, tolerance, gentle, |cloth| {
                        let at = |index: usize| cloth.positions[index];
                        segment_segment_distance(at(a), at(b), at(c), at(d))
                    });
                    stats.ccd_hits += 1;
                    clean = false;
                }
            }
            if clean {
                settled = true;
                break;
            }
        }
        self.proximity = near;
        if !settled {
            // Nicht zur Ruhe gekommen: alles zurücksetzen.
            let all: Vec<usize> = (0..self.positions.len()).collect();
            self.roll_back(&all);
            stats.full_rollbacks += 1;
        }
    }

    /// Zählt die Stellen, an denen eine Kante ein fremdes Dreieck
    /// durchstößt.
    ///
    /// Null heißt: Das Netz durchdringt weder sich selbst noch ein
    /// Hindernis. Gedacht zum Nachprüfen, nicht für jeden Zeitschritt.
    pub fn intersections(&self) -> usize {
        let mut hash = SpatialHash::build(
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
