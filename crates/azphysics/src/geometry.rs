//! Geometrie für Kollisionen: nächste Punkte, Abstände, Schnitt-Tests.
//!
//! Die Algorithmen für nächste Punkte folgen C. Ericson, „Real-Time
//! Collision Detection“ (2005), Abschnitte 5.1.2, 5.1.5 und 5.1.9.

use azmath::Vector;
use azmath::scalar::RealElementary;

/// Punkt oder Richtung im Raum.
pub type Vec3<T> = Vector<T, 3>;

/// Achsenparalleler Hüllquader in `D` Dimensionen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Aabb<T, const D: usize> {
    /// Ecke mit den kleinsten Koordinaten.
    pub min: Vector<T, D>,
    /// Ecke mit den größten Koordinaten.
    pub max: Vector<T, D>,
}

impl<T: RealElementary, const D: usize> Aabb<T, D> {
    /// Kleinster Quader, der alle Punkte enthält; `None` ohne Punkte.
    pub fn from_points(points: impl IntoIterator<Item = Vector<T, D>>) -> Option<Self> {
        let mut points = points.into_iter();
        let first = points.next()?;
        Some(points.fold(Self { min: first, max: first }, |hull, p| Self {
            min: hull.min.zip_map(p, |a, b| a.min(b)),
            max: hull.max.zip_map(p, |a, b| a.max(b)),
        }))
    }

    /// Nach allen Seiten um `margin` vergrößert.
    pub fn expanded(self, margin: T) -> Self {
        Self { min: self.min.map(|a| a - margin), max: self.max.map(|a| a + margin) }
    }

    /// Kleinster Quader, der beide enthält.
    pub fn union(self, other: Self) -> Self {
        Self {
            min: self.min.zip_map(other.min, |a, b| a.min(b)),
            max: self.max.zip_map(other.max, |a, b| a.max(b)),
        }
    }

    /// `true`, wenn sich die Quader berühren oder überlappen.
    pub fn overlaps(self, other: Self) -> bool {
        self.min.iter().zip(other.max.iter()).all(|(&lo, &hi)| lo <= hi)
            && other.min.iter().zip(self.max.iter()).all(|(&lo, &hi)| lo <= hi)
    }

    /// `true`, wenn der Punkt im Quader liegt.
    pub fn contains(self, point: Vector<T, D>) -> bool {
        self.min.iter().zip(point.iter()).all(|(&lo, &x)| lo <= x)
            && point.iter().zip(self.max.iter()).all(|(&x, &hi)| x <= hi)
    }
}

/// Auf [0, 1] begrenzt.
#[inline]
fn clamp01<T: RealElementary>(value: T) -> T {
    value.max(T::ZERO).min(T::ONE)
}

/// Nächster Punkt der Strecke `a`–`b` zu `p`.
///
/// Liefert den Punkt und den Parameter `t` in [0, 1] mit
/// `Punkt = a + t·(b − a)`.
pub fn closest_point_on_segment<T: RealElementary, const D: usize>(
    p: Vector<T, D>,
    a: Vector<T, D>,
    b: Vector<T, D>,
) -> (Vector<T, D>, T) {
    let ab = b - a;
    let length_sqr = ab.dot(ab);
    let t = if length_sqr > T::ZERO { clamp01((p - a).dot(ab) / length_sqr) } else { T::ZERO };
    (a + ab * t, t)
}

/// Nächster Punkt des Dreiecks `a`, `b`, `c` zu `p`.
///
/// Liefert den Punkt und seine baryzentrischen Koordinaten `[u, v, w]` mit
/// `Punkt = u·a + v·b + w·c`.
pub fn closest_point_on_triangle<T: RealElementary>(
    p: Vec3<T>,
    a: Vec3<T>,
    b: Vec3<T>,
    c: Vec3<T>,
) -> (Vec3<T>, [T; 3]) {
    let (zero, one) = (T::ZERO, T::ONE);
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= zero && d2 <= zero {
        return (a, [one, zero, zero]);
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= zero && d4 <= d3 {
        return (b, [zero, one, zero]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= zero && d1 >= zero && d3 <= zero && d1 - d3 > zero {
        let v = d1 / (d1 - d3);
        return (a + ab * v, [one - v, v, zero]);
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= zero && d5 <= d6 {
        return (c, [zero, zero, one]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= zero && d2 >= zero && d6 <= zero && d2 - d6 > zero {
        let w = d2 / (d2 - d6);
        return (a + ac * w, [one - w, zero, w]);
    }
    let va = d3 * d6 - d5 * d4;
    let (e1, e2) = (d4 - d3, d5 - d6);
    if va <= zero && e1 >= zero && e2 >= zero && e1 + e2 > zero {
        let w = e1 / (e1 + e2);
        return (b + (c - b) * w, [zero, one - w, w]);
    }
    let sum = va + vb + vc;
    if sum > zero {
        let (v, w) = (vb / sum, vc / sum);
        return (a + ab * v + ac * w, [one - v - w, v, w]);
    }
    // Entartetes Dreieck (Fläche null): die nächste der drei Kanten.
    let candidates = [
        (closest_point_on_segment(p, a, b), 0),
        (closest_point_on_segment(p, b, c), 1),
        (closest_point_on_segment(p, c, a), 2),
    ];
    let mut best = candidates[0];
    for candidate in candidates {
        if (p - candidate.0.0).norm_sqr() < (p - best.0.0).norm_sqr() {
            best = candidate;
        }
    }
    let ((point, t), edge) = best;
    let weights = match edge {
        0 => [one - t, t, zero],
        1 => [zero, one - t, t],
        _ => [t, zero, one - t],
    };
    (point, weights)
}

/// Nächste Punkte der Strecken `p1`–`q1` und `p2`–`q2`.
///
/// Liefert beide Punkte und ihre Parameter `s`, `t` in [0, 1].
pub fn closest_points_on_segments<T: RealElementary>(
    p1: Vec3<T>,
    q1: Vec3<T>,
    p2: Vec3<T>,
    q2: Vec3<T>,
) -> (Vec3<T>, Vec3<T>, T, T) {
    let (d1, d2, r) = (q1 - p1, q2 - p2, p1 - p2);
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let zero = T::ZERO;
    let (s, t) = if a <= zero && e <= zero {
        (zero, zero)
    } else if a <= zero {
        (zero, clamp01(f / e))
    } else {
        let c = d1.dot(r);
        if e <= zero {
            (clamp01(-c / a), zero)
        } else {
            let b = d1.dot(d2);
            let denominator = a * e - b * b;
            let s = if denominator > zero { clamp01((b * f - c * e) / denominator) } else { zero };
            let t = (b * s + f) / e;
            if t < zero {
                (clamp01(-c / a), zero)
            } else if t > T::ONE {
                (clamp01((b - c) / a), T::ONE)
            } else {
                (s, t)
            }
        }
    };
    (p1 + d1 * s, p2 + d2 * t, s, t)
}

/// Abstand des Punkts `p` vom Dreieck `a`, `b`, `c`.
pub fn point_triangle_distance<T: RealElementary>(p: Vec3<T>, a: Vec3<T>, b: Vec3<T>, c: Vec3<T>) -> T {
    (p - closest_point_on_triangle(p, a, b, c).0).norm()
}

/// Abstand der Strecken `p1`–`q1` und `p2`–`q2`.
pub fn segment_segment_distance<T: RealElementary>(
    p1: Vec3<T>,
    q1: Vec3<T>,
    p2: Vec3<T>,
    q2: Vec3<T>,
) -> T {
    let (c1, c2, _, _) = closest_points_on_segments(p1, q1, p2, q2);
    (c1 - c2).norm()
}

/// Nicht normierte Normale des Dreiecks: `(b − a) × (c − a)`; ihre Länge ist
/// die doppelte Fläche.
#[inline]
pub fn triangle_normal<T: RealElementary>(a: Vec3<T>, b: Vec3<T>, c: Vec3<T>) -> Vec3<T> {
    (b - a).cross(c - a)
}

/// `true`, wenn die Strecke `p`–`q` das Dreieck `a`, `b`, `c` schneidet.
///
/// Berührungen genau am Rand oder Strecken in der Dreiecksebene gelten
/// nicht als Schnitt. Gedacht zum Nachprüfen, ob ein Netz sich selbst
/// durchdringt.
pub fn segment_intersects_triangle<T: RealElementary>(
    p: Vec3<T>,
    q: Vec3<T>,
    a: Vec3<T>,
    b: Vec3<T>,
    c: Vec3<T>,
) -> bool {
    let normal = triangle_normal(a, b, c);
    let (side_p, side_q) = (normal.dot(p - a), normal.dot(q - a));
    // Beide Enden auf derselben Seite oder in der Ebene: kein Schnitt.
    if side_p * side_q >= T::ZERO {
        return false;
    }
    let t = side_p / (side_p - side_q);
    let hit = p + (q - p) * t;
    // Liegt der Durchstoßpunkt echt im Inneren?
    let inside = |u: Vec3<T>, v: Vec3<T>| (v - u).cross(hit - u).dot(normal) > T::ZERO;
    inside(a, b) && inside(b, c) && inside(c, a)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vec3<f64> {
        Vector([x, y, z])
    }

    fn near(a: Vec3<f64>, b: Vec3<f64>) -> bool {
        (a - b).norm() < 1.0e-12
    }

    #[test]
    fn triangle_regions() {
        let (a, b, c) = (v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0));
        // Über dem Inneren: senkrechte Projektion.
        let (point, weights) = closest_point_on_triangle(v(0.5, 0.5, 3.0), a, b, c);
        assert!(near(point, v(0.5, 0.5, 0.0)));
        assert!(near(a * weights[0] + b * weights[1] + c * weights[2], point));
        assert!((point_triangle_distance(v(0.5, 0.5, 3.0), a, b, c) - 3.0).abs() < 1.0e-12);
        // Ecken, Kanten.
        assert!(near(closest_point_on_triangle(v(-1.0, -1.0, 0.0), a, b, c).0, a));
        assert!(near(closest_point_on_triangle(v(3.0, -1.0, 1.0), a, b, c).0, b));
        assert!(near(closest_point_on_triangle(v(-1.0, 3.0, 0.0), a, b, c).0, c));
        assert!(near(closest_point_on_triangle(v(1.0, -1.0, 0.0), a, b, c).0, v(1.0, 0.0, 0.0)));
        assert!(near(closest_point_on_triangle(v(-1.0, 1.0, 0.0), a, b, c).0, v(0.0, 1.0, 0.0)));
        assert!(near(closest_point_on_triangle(v(2.0, 2.0, 0.0), a, b, c).0, v(1.0, 1.0, 0.0)));
        // Baryzentrische Koordinaten summieren sich immer zu eins.
        for p in [v(0.3, 0.2, 1.0), v(5.0, 5.0, 5.0), v(-2.0, 0.5, -1.0), v(1.0, -3.0, 0.0)] {
            let (point, w) = closest_point_on_triangle(p, a, b, c);
            assert!((w[0] + w[1] + w[2] - 1.0).abs() < 1.0e-12);
            assert!(w.iter().all(|&weight| weight >= -1.0e-12));
            assert!(near(a * w[0] + b * w[1] + c * w[2], point));
        }
        // Entartet: alle Ecken auf einer Geraden.
        let (point, _) = closest_point_on_triangle(v(1.0, 1.0, 0.0), a, b, v(4.0, 0.0, 0.0));
        assert!(near(point, v(1.0, 0.0, 0.0)));
    }

    #[test]
    fn segments() {
        // Windschief gekreuzt: Abstand ist der Höhenunterschied.
        let (c1, c2, s, t) = closest_points_on_segments(
            v(-1.0, 0.0, 0.0),
            v(1.0, 0.0, 0.0),
            v(0.0, -1.0, 2.0),
            v(0.0, 1.0, 2.0),
        );
        assert!(near(c1, v(0.0, 0.0, 0.0)) && near(c2, v(0.0, 0.0, 2.0)));
        assert!((s - 0.5).abs() < 1.0e-12 && (t - 0.5).abs() < 1.0e-12);
        // Parallel, Ende an Ende, zu Punkten entartet.
        let d = segment_segment_distance(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(1.0, 1.0, 0.0));
        assert!((d - 1.0).abs() < 1.0e-12);
        let d = segment_segment_distance(v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(3.0, 0.0, 0.0), v(5.0, 0.0, 0.0));
        assert!((d - 2.0).abs() < 1.0e-12);
        let d = segment_segment_distance(v(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0), v(3.0, 4.0, 0.0), v(3.0, 4.0, 0.0));
        assert!((d - 5.0).abs() < 1.0e-12);
        let (point, t) = closest_point_on_segment(v(0.25, 1.0, 0.0), v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0));
        assert!(near(point, v(0.25, 0.0, 0.0)) && (t - 0.25).abs() < 1.0e-12);
        let (point, t) = closest_point_on_segment(v(5.0, 1.0, 0.0), v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0));
        assert!(near(point, v(1.0, 0.0, 0.0)) && t == 1.0);
    }

    #[test]
    fn intersections_and_boxes() {
        let (a, b, c) = (v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0));
        assert!(segment_intersects_triangle(v(0.5, 0.5, -1.0), v(0.5, 0.5, 1.0), a, b, c));
        assert!(segment_intersects_triangle(v(0.5, 0.5, 1.0), v(0.5, 0.5, -1.0), a, b, c));
        assert!(!segment_intersects_triangle(v(0.5, 0.5, 1.0), v(0.5, 0.5, 2.0), a, b, c));
        assert!(!segment_intersects_triangle(v(3.0, 3.0, -1.0), v(3.0, 3.0, 1.0), a, b, c));
        assert!(!segment_intersects_triangle(v(0.0, 0.0, 0.0), v(1.0, 1.0, 0.0), a, b, c));
        assert!((triangle_normal(a, b, c) - v(0.0, 0.0, 4.0)).norm() < 1.0e-12);

        let hull = Aabb::from_points([a, b, c]);
        assert!(hull.is_some_and(|h| {
            h.min == v(0.0, 0.0, 0.0)
                && h.max == v(2.0, 2.0, 0.0)
                && h.contains(v(1.0, 1.0, 0.0))
                && !h.contains(v(1.0, 1.0, 0.1))
                && h.expanded(0.5).contains(v(1.0, 1.0, 0.1))
                && h.overlaps(Aabb { min: v(2.0, 2.0, 0.0), max: v(3.0, 3.0, 1.0) })
                && !h.overlaps(Aabb { min: v(2.1, 0.0, 0.0), max: v(3.0, 3.0, 1.0) })
                && h.union(Aabb { min: v(-1.0, 0.0, 0.0), max: v(0.0, 5.0, 1.0) }).max == v(2.0, 5.0, 1.0)
        }));
        assert!(Aabb::<f64, 3>::from_points([]).is_none());
    }
}
