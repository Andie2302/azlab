//! Kollision des Tuchs mit Dreiecken: Abstoßung und Sicherung gegen
//! Durchdringung.

use super::model::{Cloth, StepStats, StepPhase};
use super::par::gather;
use crate::broadphase::{QueryScratch, SpatialHash};
use crate::ccd::{Contact, edge_edge, vertex_triangle};
use crate::geometry::{
    Aabb, Vec3, closest_point_on_triangle, point_triangle_distance, segment_intersects_triangle,
    segment_segment_distance, triangle_normal,
};
use alloc::vec::Vec;
use azmath::scalar::RealElementary;

/// Durchläufe der Sicherung, bevor alles zurückgesetzt wird.
const MAX_PASSES: usize = 8;

/// Sicherheitsabstand der Nachbarschaftslisten in Maschenweiten. Größer
/// heißt seltener neu bauen, aber mehr Paare in jedem Teilschritt.
const SLACK: f64 = 0.25;

impl<T: RealElementary + Send + Sync> Cloth<T> {
    /// Kantenlänge der Suchzellen: etwa eine Maschenweite.
    fn cell_size(&self) -> T {
        let (sum, count) = self.links.iter().filter(|link| !link.bend).fold(
            (T::ZERO, T::ZERO),
            |(sum, count), link| (sum + link.rest, count + T::ONE),
        );
        let mesh = if count > T::ZERO { sum / count } else { T::ONE };
        mesh.max(self.params.thickness * T::from_f64(2.0))
    }

    /// Abstand, den die Sicherung als Berührung wertet.
    fn ccd_tolerance(&self) -> T {
        self.params.thickness * T::from_f64(0.1)
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
    fn refresh_proximity(&mut self, stats: &mut StepStats) {
        let thickness = self.params.thickness;
        let near = &self.proximity;
        let mut current =
            near.valid && near.thickness == thickness && near.anchors.len() == self.positions.len();
        if current {
            // Bewegt sich das Tuch als Ganzes, ändert das nichts daran,
            // welche Teile sich nahe sind. Die gemeinsame Verschiebung seit
            // dem Bau wird deshalb herausgerechnet – sonst würden die Listen
            // im freien Fall in fast jedem Teilschritt neu gebaut. Das geht
            // nur, wenn sich alle Teilchen bewegen können: Gegenüber einem
            // festen Teilchen ist die Verschiebung eine echte Annäherung.
            let shift = if near.all_movable && !self.positions.is_empty() {
                let sum = near.anchors.iter().zip(&self.positions).fold(Vec3::ZERO, |sum, (anchor, &position)| {
                    sum + (position - (anchor.min + anchor.max) * T::from_f64(0.5))
                });
                sum * T::from_f64(self.positions.len() as f64).recip()
            } else {
                Vec3::ZERO
            };
            current = near
                .anchors
                .iter()
                .zip(&self.positions)
                .all(|(anchor, &position)| anchor.expanded(near.slack).contains(position - shift));
        }
        if current {
            return;
        }

        let mesh = self.cell_size();
        let slack = mesh * T::from_f64(SLACK);
        let mut near = core::mem::take(&mut self.proximity);
        near.anchors.clear();
        near.anchors.extend(
            self.previous.iter().zip(&self.positions).map(|(&from, &to)| Aabb {
                min: from.zip_map(to, |a, b| a.min(b)),
                max: from.zip_map(to, |a, b| a.max(b)),
            }),
        );
        // Bis zum nächsten Bau bleibt jedes Teilchen in seinem Bereich,
        // vergrößert um `slack`. Ein Paar kann sich also nur dann näher als
        // `contact` kommen, wenn die Bereiche des einen, vergrößert um
        // `slack`, und des anderen, vergrößert um `slack + contact`,
        // einander berühren.
        let hull = |indices: &[usize], margin: T| {
            indices[1..].iter().fold(near.anchors[indices[0]], |hull, &i| hull.union(near.anchors[i])).expanded(margin)
        };
        let two = T::from_f64(2.0);

        // Teilchen gegen Dreieck: gebraucht von der Abstoßung (Abstand bis
        // zur Dicke) und von der Sicherung (Abstand bis zur Toleranz).
        let reach = slack + thickness;
        let triangle_boxes: Vec<Aabb<T, 3>> = self.triangles.iter().map(|tri| hull(tri, reach)).collect();
        // Zellen so groß wie die Quader, die hineinsortiert werden: Jeder
        // fällt dann in wenige Zellen. Kleinere Zellen vervielfachen die
        // Einträge und machen den Bau um ein Vielfaches teurer.
        let hash = SpatialHash::build(mesh + reach * two, triangle_boxes.iter().copied().enumerate());
        // Die Suche je Teilchen ist unabhängig von der der anderen und
        // lässt sich auf mehrere Kerne verteilen.
        near.vertex_triangle = gather(self.positions.len(), |vertices, found| {
            let (mut scratch, mut candidates) = (QueryScratch::default(), Vec::new());
            for vertex in vertices {
                let bounds = hull(&[vertex], slack);
                hash.query_with(&mut scratch, bounds, &mut candidates);
                for &id in &candidates {
                    let triangle = self.triangles[id];
                    if !triangle.contains(&vertex)
                        && triangle_boxes[id].overlaps(bounds)
                        && !(self.all_fixed(&[vertex]) && self.all_fixed(&triangle))
                    {
                        found.push((vertex, id));
                    }
                }
            }
        });

        // Kante gegen Kante: nur von der Sicherung gebraucht, deshalb mit
        // der kleineren Reichweite. Beide Kanten bekommen die halbe
        // Toleranz dazu; so ist die Prüfung in beiden Richtungen dieselbe
        // und jedes Paar wird nur einmal betrachtet.
        let reach = slack + self.ccd_tolerance() * T::from_f64(0.5);
        let edge_boxes: Vec<Aabb<T, 3>> = self.edges.iter().map(|edge| hull(edge, reach)).collect();
        let hash = SpatialHash::build(mesh + reach * two, edge_boxes.iter().copied().enumerate());
        near.edge_edge = gather(self.edges.len(), |edges, found| {
            let (mut scratch, mut candidates) = (QueryScratch::default(), Vec::new());
            for first in edges {
                let [a, b] = self.edges[first];
                hash.query_with(&mut scratch, edge_boxes[first], &mut candidates);
                for &second in &candidates {
                    let [c, d] = self.edges[second];
                    // Kanten mit gemeinsamem Teilchen können sich nicht
                    // durchdringen.
                    if second > first
                        && a != c
                        && a != d
                        && b != c
                        && b != d
                        && edge_boxes[first].overlaps(edge_boxes[second])
                        && !self.all_fixed(&[a, b, c, d])
                    {
                        found.push((first, second));
                    }
                }
            }
        });
        stats.proximity_rebuilds += 1;
        near.slack = slack;
        near.thickness = thickness;
        near.all_movable = self.inverse_mass.iter().all(|&inverse| inverse > T::ZERO);
        near.valid = true;
        self.proximity = near;
    }

    /// `true`, wenn das Teilchen dem Dreieck näher ist als `distance`.
    fn closer_than(&self, vertex: usize, id: usize, distance: T) -> bool {
        let [a, b, c] = self.triangles[id].map(|index| self.positions[index]);
        let point = self.positions[vertex];
        // Schnelle Vorprüfung: außerhalb des Hüllquaders ist das Teilchen
        // weiter weg.
        let inside = (0..3).all(|axis| {
            let (low, high) = (a[axis].min(b[axis]).min(c[axis]), a[axis].max(b[axis]).max(c[axis]));
            point[axis] >= low - distance && point[axis] <= high + distance
        });
        inside && (point - closest_point_on_triangle(point, a, b, c).0).norm() < distance
    }

    /// Schiebt Teilchen und Dreiecke auseinander, die sich näher kommen als
    /// die Dicke des Tuchs.
    pub(super) fn repel(&mut self, stats: &mut StepStats, observe: &mut impl FnMut(StepPhase)) {
        let thickness = self.params.thickness;
        self.refresh_proximity(stats);
        observe(StepPhase::Proximity);
        let near = core::mem::take(&mut self.proximity);
        // Erkennen: welche Paare sind sich zu Beginn näher als die Dicke?
        // Das ändert nichts und lässt sich auf mehrere Kerne verteilen;
        // meist sind es wenige.
        let close = gather(near.vertex_triangle.len(), |pairs, found| {
            for pair in pairs {
                let (vertex, id) = near.vertex_triangle[pair];
                if self.closer_than(vertex, id, thickness) {
                    found.push(pair);
                }
            }
        });
        // Auseinanderschieben: der Reihe nach, jedes Paar gegen den
        // aktuellen Stand. Was dabei erst zu nahe gerät, ist im nächsten
        // Teilschritt an der Reihe.
        for pair in close {
            let (vertex, id) = near.vertex_triangle[pair];
            let [a, b, c] = self.triangles[id];
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
            // Reibung nach Coulomb, wie an Ebene und Kugel: Die Strecke, um
            // die auseinandergeschoben wurde, steht für die Kraft, mit der
            // Teilchen und Dreieck aufeinander drücken. Das Gleiten des
            // Teilchens gegenüber seinem Berührpunkt auf dem Dreieck wird um
            // höchstens Reibwert mal diese Strecke gekürzt; gleitet es
            // weniger weit, haftet es.
            let friction = self.params.friction.max(T::ZERO);
            if friction > T::ZERO {
                let travel = |index: usize| self.positions[index] - self.previous[index];
                let relative =
                    travel(vertex) - travel(a) * weights[0] - travel(b) * weights[1] - travel(c) * weights[2];
                let along = relative - direction * relative.dot(direction);
                let slip = along.norm();
                if slip > T::ZERO {
                    let cut = along * (slip.min(friction * gap) / slip / share);
                    self.positions[vertex] -= cut * inverse[0];
                    self.positions[a] += cut * (weights[0] * inverse[1]);
                    self.positions[b] += cut * (weights[1] * inverse[2]);
                    self.positions[c] += cut * (weights[2] * inverse[3]);
                }
            }
            stats.repulsions += 1;
        }
        self.proximity = near;
    }

    /// Setzt Teilchen auf ihre Lage vom Anfang des Teilschritts zurück.
    fn roll_back(&mut self, indices: &[usize]) {
        for &index in indices {
            self.positions[index] = self.previous[index];
            self.carry[index] = self.previous_carry[index];
            self.halted[index] = true;
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

    /// Prüft die Bewegung eines Teilchens gegen die eines Dreiecks.
    fn vertex_contact(&self, vertex: usize, id: usize, tolerance: T) -> Contact<T> {
        let [a, b, c] = self.triangles[id];
        if [vertex, a, b, c].iter().all(|&i| self.positions[i] == self.previous[i]) {
            // Nichts bewegt sich: keine neue Berührung.
            return Contact::Free;
        }
        vertex_triangle(self.motion(vertex), self.motion(a), self.motion(b), self.motion(c), tolerance)
    }

    /// Prüft die Bewegungen zweier Kanten gegeneinander.
    fn edge_contact(&self, first: usize, second: usize, tolerance: T) -> Contact<T> {
        let ([a, b], [c, d]) = (self.edges[first], self.edges[second]);
        if [a, b, c, d].iter().all(|&i| self.positions[i] == self.previous[i]) {
            return Contact::Free;
        }
        edge_edge(self.motion(a), self.motion(b), self.motion(c), self.motion(d), tolerance)
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
        let tolerance = self.ccd_tolerance();
        self.refresh_proximity(stats);
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
            // Erkennen: alle Paare gegen den Stand zu Beginn des Durchlaufs
            // prüfen. Das ändert nichts und lässt sich auf mehrere Kerne
            // verteilen; meist findet es nichts.
            let touching_triangles = gather(near.vertex_triangle.len(), |pairs, found| {
                for pair in pairs {
                    let (vertex, id) = near.vertex_triangle[pair];
                    if triangle_boxes[id].overlaps(paths[vertex]) && !self.vertex_contact(vertex, id, tolerance).is_free() {
                        found.push(pair);
                    }
                }
            });
            let touching_edges = gather(near.edge_edge.len(), |pairs, found| {
                for pair in pairs {
                    let (first, second) = near.edge_edge[pair];
                    if edge_boxes[first].overlaps(edge_boxes[second])
                        && !self.edge_contact(first, second, tolerance).is_free()
                    {
                        found.push(pair);
                    }
                }
            });
            if touching_triangles.is_empty() && touching_edges.is_empty() {
                settled = true;
                break;
            }
            // Beheben: der Reihe nach. Jede Kürzung ändert die Bewegung der
            // Nachbarn, deshalb wird jedes gemeldete Paar hier noch einmal
            // gegen den aktuellen Stand geprüft. Was eine Kürzung an neuen
            // Berührungen auslöst, findet der nächste Durchlauf.
            for pair in touching_triangles {
                let (vertex, id) = near.vertex_triangle[pair];
                let hit = self.vertex_contact(vertex, id, tolerance);
                // Unentschieden zählt wie eine Berührung: Die Bewegung endet
                // vor dem Zeitpunkt, bis zu dem sie nachweislich frei ist.
                if let Some(contact) = hit.safe_until() {
                    let [a, b, c] = self.triangles[id];
                    self.retreat(&[vertex, a, b, c], contact, tolerance, gentle, |cloth| {
                        let at = |index: usize| cloth.positions[index];
                        point_triangle_distance(at(vertex), at(a), at(b), at(c))
                    });
                    stats.ccd_hits += 1;
                    stats.ccd_undecided += u32::from(matches!(hit, Contact::Undecided(_)));
                }
            }
            for pair in touching_edges {
                let (first, second) = near.edge_edge[pair];
                let hit = self.edge_contact(first, second, tolerance);
                if let Some(contact) = hit.safe_until() {
                    let ([a, b], [c, d]) = (self.edges[first], self.edges[second]);
                    self.retreat(&[a, b, c, d], contact, tolerance, gentle, |cloth| {
                        let at = |index: usize| cloth.positions[index];
                        segment_segment_distance(at(a), at(b), at(c), at(d))
                    });
                    stats.ccd_hits += 1;
                    stats.ccd_undecided += u32::from(matches!(hit, Contact::Undecided(_)));
                }
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
