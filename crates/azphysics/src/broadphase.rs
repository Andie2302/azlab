//! Suchstruktur für die Grobphase der Kollisionserkennung.
//!
//! Jedes Paar aus Punkt und Dreieck einzeln zu prüfen, kostet quadratisch
//! viel. Ein gleichmäßiges Raumgitter sortiert die Objekte nach ihrer Lage
//! vor: Nur was in dieselben Zellen fällt, kann sich berühren.

use crate::geometry::Aabb;
use alloc::vec::Vec;
use azmath::scalar::RealElementary;

/// Objekte, die mehr Zellen überdecken, werden gesondert geführt und bei
/// jeder Anfrage gemeldet.
const MAX_CELLS: i64 = 512;

/// Gleichmäßiges Raumgitter über Hüllquadern.
///
/// Die Suche ist vorsichtig: Sie meldet jedes Objekt, dessen Hüllquader den
/// angefragten Quader berührt, und gelegentlich zusätzliche. Sie lässt
/// keines aus.
#[derive(Clone, Debug)]
pub struct SpatialHash<T> {
    cell: T,
    /// Streutabelle über die Zellen: Die Objekte des Fachs `b` stehen in
    /// `items[starts[b]..starts[b + 1]]`.
    starts: Vec<usize>,
    items: Vec<usize>,
    /// Zu jedem Eintrag in `items` der Schlüssel seiner Zelle.
    keys: Vec<u64>,
    large: Vec<usize>,
    /// Je Kennnummer die Anfrage, in der sie zuletzt gemeldet wurde.
    seen: Vec<u32>,
    query: u32,
}

impl<T: RealElementary> SpatialHash<T> {
    /// Baut das Gitter mit Zellen der Kantenlänge `cell` aus Paaren von
    /// Kennnummer und Hüllquader.
    pub fn build(cell: T, items: impl IntoIterator<Item = (usize, Aabb<T, 3>)>) -> Self {
        let cell = if cell > T::ZERO && cell.is_finite() { cell } else { T::ONE };
        let mut hash =
            Self {
            cell,
            starts: Vec::new(),
            items: Vec::new(),
            keys: Vec::new(),
            large: Vec::new(),
            seen: Vec::new(),
            query: 0,
        };
        // Erst alle Einträge (Zellschlüssel, Kennnummer) sammeln, dann nach
        // Fächern ordnen: Zählsortierung, ohne Vergleiche.
        let mut entries: Vec<(u64, usize)> = Vec::new();
        let mut ids = 0;
        for (id, bounds) in items {
            ids = ids.max(id + 1);
            match hash.cells(bounds) {
                Some((lo, hi)) => {
                    for x in lo[0]..=hi[0] {
                        for y in lo[1]..=hi[1] {
                            for z in lo[2]..=hi[2] {
                                entries.push((key(x, y, z), id));
                            }
                        }
                    }
                }
                None => hash.large.push(id),
            }
        }
        let buckets = (entries.len() * 2).next_power_of_two().max(16);
        hash.starts.resize(buckets + 1, 0);
        for &(key, _) in &entries {
            hash.starts[(key as usize & (buckets - 1)) + 1] += 1;
        }
        for bucket in 0..buckets {
            hash.starts[bucket + 1] += hash.starts[bucket];
        }
        let mut next = hash.starts.clone();
        hash.items.resize(entries.len(), 0);
        hash.keys.resize(entries.len(), 0);
        for &(key, id) in &entries {
            let slot = &mut next[key as usize & (buckets - 1)];
            hash.items[*slot] = id;
            hash.keys[*slot] = key;
            *slot += 1;
        }
        hash.seen.resize(ids, 0);
        hash
    }

    /// Zellbereich eines Quaders; `None`, wenn er zu viele Zellen überdeckt
    /// oder nicht endlich ist.
    fn cells(&self, bounds: Aabb<T, 3>) -> Option<([i64; 3], [i64; 3])> {
        let index = |value: T| {
            let scaled = (value / self.cell).floor();
            scaled.is_finite().then(|| scaled.to_f64() as i64)
        };
        let mut lo = [0_i64; 3];
        let mut hi = [0_i64; 3];
        let mut count = 1_i64;
        for axis in 0..3 {
            lo[axis] = index(bounds.min[axis])?;
            hi[axis] = index(bounds.max[axis])?;
            count = count.saturating_mul(hi[axis].saturating_sub(lo[axis]).saturating_add(1));
        }
        (count > 0 && count <= MAX_CELLS).then_some((lo, hi))
    }

    /// Schreibt die Kennnummern aller Objekte nach `out`, deren Hüllquader
    /// `bounds` berühren könnte; jede einmal, in keiner bestimmten
    /// Reihenfolge.
    pub fn query(&mut self, bounds: Aabb<T, 3>, out: &mut Vec<usize>) {
        out.clear();
        out.extend_from_slice(&self.large);
        // Neue Anfragenummer; beim Überlauf die Merkliste leeren.
        self.query = self.query.wrapping_add(1);
        if self.query == 0 {
            self.seen.fill(0);
            self.query = 1;
        }
        let buckets = self.starts.len() - 1;
        let range = match self.cells(bounds) {
            Some(range) => range,
            // Riesiger Suchbereich: alles melden.
            None => {
                for &id in &self.items {
                    if self.seen[id] != self.query {
                        self.seen[id] = self.query;
                        out.push(id);
                    }
                }
                return;
            }
        };
        let (lo, hi) = range;
        for x in lo[0]..=hi[0] {
            for y in lo[1]..=hi[1] {
                for z in lo[2]..=hi[2] {
                    let wanted = key(x, y, z);
                    let bucket = wanted as usize & (buckets - 1);
                    for slot in self.starts[bucket]..self.starts[bucket + 1] {
                        let id = self.items[slot];
                        if self.keys[slot] == wanted && self.seen[id] != self.query {
                            self.seen[id] = self.query;
                            out.push(id);
                        }
                    }
                }
            }
        }
    }
}

/// Schlüssel einer Zelle. Verschiedene Zellen können denselben Schlüssel
/// bekommen; das führt nur zu überzähligen, nie zu fehlenden Treffern.
#[inline]
fn key(x: i64, y: i64, z: i64) -> u64 {
    (x as u64).wrapping_mul(73_856_093)
        ^ (y as u64).wrapping_mul(19_349_663)
        ^ (z as u64).wrapping_mul(83_492_791)
}

#[cfg(all(test, any(feature = "std", feature = "libm")))]
mod tests {
    use super::*;
    use alloc::vec;
    use azmath::Vector;

    fn cube(center: [f64; 3], half: f64) -> Aabb<f64, 3> {
        Aabb { min: Vector(center.map(|c| c - half)), max: Vector(center.map(|c| c + half)) }
    }

    #[test]
    fn finds_every_overlap() {
        // Würfel auf einem Raster; jede Anfrage muss alle überlappenden finden.
        let boxes: Vec<Aabb<f64, 3>> = (0..125)
            .map(|i| cube([(i % 5) as f64, ((i / 5) % 5) as f64, (i / 25) as f64], 0.3))
            .collect();
        let mut hash = SpatialHash::build(0.75, boxes.iter().copied().enumerate());
        let mut found = Vec::new();
        for probe in [cube([2.0, 2.0, 2.0], 0.5), cube([0.1, 3.9, 1.5], 1.2), cube([9.0, 9.0, 9.0], 0.2)] {
            hash.query(probe, &mut found);
            for (id, bounds) in boxes.iter().enumerate() {
                if bounds.overlaps(probe) {
                    assert!(found.contains(&id));
                }
            }
            // Keine Kennnummer doppelt.
            let mut sorted = found.clone();
            sorted.sort_unstable();
            assert!(sorted.windows(2).all(|pair| pair[0] < pair[1]));
        }
        // Weit weg: nichts (bis auf mögliche Zufallstreffer, hier keine).
        hash.query(cube([100.0, 100.0, 100.0], 0.1), &mut found);
        assert!(found.len() < 5);
    }

    #[test]
    fn large_and_degenerate() {
        let huge = cube([0.0, 0.0, 0.0], 1000.0);
        let small = cube([5.0, 5.0, 5.0], 0.1);
        let mut hash = SpatialHash::build(0.5, [(7, huge), (3, small)]);
        let mut found = vec![99];
        // Das riesige Objekt wird immer gemeldet.
        hash.query(cube([50.0, 0.0, 0.0], 0.1), &mut found);
        assert_eq!(found, vec![7]);
        hash.query(small, &mut found);
        found.sort_unstable();
        assert_eq!(found, vec![3, 7]);
        // Riesige Anfrage meldet alles.
        hash.query(huge, &mut found);
        found.sort_unstable();
        assert_eq!(found, vec![3, 7]);
        // Ungültige Zellgröße wird ersetzt, nicht endliche Quader gelten als groß.
        let mut odd = SpatialHash::build(f64::NAN, [(1, Aabb { min: Vector([f64::NAN; 3]), max: Vector([1.0; 3]) })]);
        odd.query(small, &mut found);
        assert_eq!(found, vec![1]);
    }
}
