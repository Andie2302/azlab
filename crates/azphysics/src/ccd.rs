//! Kontinuierliche Kollisionserkennung (CCD).
//!
//! Gewöhnliche Kollisionsprüfungen sehen nur die Lage am Ende eines
//! Zeitschritts. Bewegt sich etwas schnell oder ist es dünn, kann es
//! zwischen zwei Schritten durch ein Hindernis hindurchtreten, ohne dass je
//! eine Überlappung sichtbar wird („Tunneln“). Die kontinuierliche Prüfung
//! betrachtet stattdessen die ganze Bewegung während des Schritts.
//!
//! Jeder Punkt bewegt sich im Schritt geradlinig von seiner Anfangs- zu
//! seiner Endlage. Für ein Dreiecksnetz gibt es genau zwei Arten, wie sich
//! zwei Teile durchdringen können:
//!
//! - ein Punkt tritt durch ein Dreieck ([`vertex_triangle`]),
//! - zwei Kanten treten durcheinander ([`edge_edge`]).
//!
//! In beiden Fällen liegen die vier beteiligten Punkte im Augenblick der
//! Berührung in einer Ebene. Diese Bedingung ist ein Polynom dritten Grades
//! in der Zeit. Seine Nullstellen im Schritt werden mit der
//! Intervall-Einschließung aus azmath gesucht, die keine Nullstelle
//! übersehen kann; an jeder wird geprüft, ob sich die Teile dort wirklich
//! treffen.
//!
//! Bleiben die vier Punkte die ganze Zeit in einer Ebene (ein flaches Tuch,
//! das in sich gleitet), sagt das Polynom nichts aus. Dann entscheidet der
//! Abstand selbst: Er wird schrittweise verfolgt, jeweils so weit, wie sich
//! die Teile in der Zeit höchstens nähern können.
//!
//! Gerechnet wird mit der Bewegung der Teile gegeneinander. Eine gemeinsame
//! Verschiebung – das Tuch fliegt als Ganzes – zählt nicht als Annäherung.
//!
//! Die Prüfung ist vorsichtig: Im Zweifel meldet sie eine Berührung. Sie
//! meldet nie „frei“, wenn sich die Teile durchdringen.

use crate::geometry::{Vec3, closest_points_on_segments, point_triangle_distance, segment_segment_distance};
use azmath::interval::Interval;
use azmath::poly::{BracketKind, RootBracket, isolate_real_roots, refine_root};
use azmath::scalar::RealElementary;
use azmath::solve::Tolerance;

/// Lage am Anfang und am Ende des Zeitschritts.
pub type Motion<T> = [Vec3<T>; 2];

/// Lage zur Zeit `t` in [0, 1].
#[inline]
fn at<T: RealElementary>(motion: &Motion<T>, t: T) -> Vec3<T> {
    motion[0] + (motion[1] - motion[0]) * t
}

/// Größte Strecke, die einer der Punkte im Schritt zurücklegt.
fn reach<T: RealElementary>(motions: &[&Motion<T>]) -> T {
    motions.iter().fold(T::ZERO, |longest, motion| longest.max((motion[1] - motion[0]).norm()))
}

/// Koeffizienten (aufsteigend) von `(e1(t) × e2(t)) · d(t)`, wenn sich
/// `e1`, `e2` und `d` geradlinig ändern.
fn coplanarity<T: RealElementary>(e1: Motion<T>, e2: Motion<T>, d: Motion<T>) -> [T; 4] {
    let (e1_0, e1_v) = (e1[0], e1[1] - e1[0]);
    let (e2_0, e2_v) = (e2[0], e2[1] - e2[0]);
    let (d_0, d_v) = (d[0], d[1] - d[0]);
    // Normale n(t) = n0 + t·n1 + t²·n2
    let n0 = e1_0.cross(e2_0);
    let n1 = e1_0.cross(e2_v) + e1_v.cross(e2_0);
    let n2 = e1_v.cross(e2_v);
    [n0.dot(d_0), n0.dot(d_v) + n1.dot(d_0), n1.dot(d_v) + n2.dot(d_0), n2.dot(d_v)]
}

/// Zieht von vier Bewegungen ihre mittlere Verschiebung ab.
///
/// Eine gemeinsame Verschiebung ändert nichts daran, ob und wann sich die
/// Teile berühren. Ohne sie sind die Wege kürzer: Der schnelle Ausschluss
/// greift öfter, und ein Tuch, das als Ganzes fliegt, stößt nicht scheinbar
/// mit sich selbst zusammen.
fn relative<T: RealElementary>(motions: [Motion<T>; 4]) -> [Motion<T>; 4] {
    let sum = motions.iter().fold(Vec3::ZERO, |sum, motion| sum + (motion[1] - motion[0]));
    let mean = sum * T::from_f64(0.25);
    motions.map(|motion| [motion[0], motion[1] - mean])
}

/// Sucht in [`from`, `to`] den frühesten Zeitpunkt, an dem `distance(t)`
/// höchstens `tolerance` beträgt, durch vorsichtiges Vorrücken.
///
/// Der Abstand ändert sich je Zeiteinheit höchstens um `reach`. Vom
/// Zeitpunkt `t` aus kann die Berührung also frühestens nach
/// `(Abstand − Toleranz) / reach` eintreten; so weit wird vorgerückt. Das
/// übersieht keine Berührung, für keine Art von Bewegung.
fn advance<T: RealElementary>(
    from: T,
    to: T,
    reach: T,
    tolerance: T,
    distance: impl Fn(T) -> T,
) -> Option<T> {
    let mut t = from;
    for _ in 0..64 {
        let gap = distance(t) - tolerance;
        if gap <= T::ZERO || gap.is_nan() {
            return Some(t);
        }
        if reach <= T::ZERO {
            return None;
        }
        t += gap / reach;
        if t > to {
            return None;
        }
    }
    // Schleichende Annäherung: vorsichtshalber als Berührung melden.
    Some(t)
}

/// Sucht den frühesten Zeitpunkt in [0, 1], an dem die Koplanaritäts-
/// bedingung erfüllt ist und `distance(t)` höchstens `tolerance` beträgt.
///
/// `reach` ist die größte Strecke, um die sich die Teile im Schritt
/// annähern können; sie begrenzt, wie stark sich der Abstand innerhalb
/// einer Zeitspanne ändern kann.
fn earliest_contact<T: RealElementary>(
    coefficients: [T; 4],
    reach: T,
    tolerance: T,
    distance: impl Fn(T) -> T,
) -> Option<T> {
    let unit = Interval::new(T::ZERO, T::ONE)?;
    let placeholder = RootBracket { interval: unit, kind: BracketKind::Unresolved };
    let mut brackets = [placeholder; 16];
    let min_width = T::EPSILON * T::from_f64(4096.0);
    let Ok(count) = isolate_real_roots(&coefficients, unit, min_width, &mut brackets) else {
        // Polynom identisch null (die vier Punkte bleiben in einer Ebene)
        // oder zu viele Kandidaten: Die Nullstellen helfen nicht weiter,
        // der Abstand selbst entscheidet.
        return advance(T::ZERO, T::ONE, reach, tolerance, distance);
    };
    let tol = Tolerance::default();
    for bracket in brackets.iter().take(count) {
        let (lo, hi) = (bracket.interval.lo().max(T::ZERO), bracket.interval.hi().min(T::ONE));
        if bracket.kind == BracketKind::Unique {
            if let Ok(root) = refine_root(&coefficients, bracket, &tol) {
                // Die Nullstelle ist bis auf Rundung bekannt.
                let slack = reach * T::EPSILON.sqrt();
                if distance(root.x) <= tolerance + slack {
                    return Some(root.x.max(T::ZERO).min(T::ONE));
                }
                continue;
            }
        }
        // Nicht aufgelöst: Die Nullstelle liegt irgendwo in der Klammer.
        if let Some(t) = advance(lo, hi, reach, tolerance, &distance) {
            return Some(t);
        }
    }
    None
}

/// Wann berührt der bewegte Punkt `p` das bewegte Dreieck `a`, `b`, `c`?
///
/// Jedes Argument ist die Lage am Anfang und am Ende des Schritts.
/// `tolerance` ist der Abstand, der als Berührung gilt (größer als null,
/// klein gegen die Abmessungen). Liefert den frühesten Zeitpunkt in [0, 1]
/// oder `None`, wenn der Punkt das Dreieck im Schritt sicher nicht berührt.
pub fn vertex_triangle<T: RealElementary>(
    p: Motion<T>,
    a: Motion<T>,
    b: Motion<T>,
    c: Motion<T>,
    tolerance: T,
) -> Option<T> {
    let [p, a, b, c] = relative([p, a, b, c]);
    // Schneller Ausschluss: Der Abstand kann höchstens um die Wege schrumpfen,
    // die Punkt und Dreieck zurücklegen.
    let travel = (p[1] - p[0]).norm() + reach(&[&a, &b, &c]);
    let start = point_triangle_distance(p[0], a[0], b[0], c[0]);
    if start > travel + tolerance {
        return None;
    }
    if start <= tolerance {
        return Some(T::ZERO);
    }
    let coefficients = coplanarity(
        [b[0] - a[0], b[1] - a[1]],
        [c[0] - a[0], c[1] - a[1]],
        [p[0] - a[0], p[1] - a[1]],
    );
    earliest_contact(coefficients, travel, tolerance, |t| {
        point_triangle_distance(at(&p, t), at(&a, t), at(&b, t), at(&c, t))
    })
}

/// Wann berühren sich die bewegten Kanten `a`–`b` und `c`–`d`?
///
/// Argumente und Ergebnis wie bei [`vertex_triangle`].
pub fn edge_edge<T: RealElementary>(
    a: Motion<T>,
    b: Motion<T>,
    c: Motion<T>,
    d: Motion<T>,
    tolerance: T,
) -> Option<T> {
    let [a, b, c, d] = relative([a, b, c, d]);
    let travel = reach(&[&a, &b]) + reach(&[&c, &d]);
    let start = segment_segment_distance(a[0], b[0], c[0], d[0]);
    if start > travel + tolerance {
        return None;
    }
    if start <= tolerance {
        return Some(T::ZERO);
    }
    let coefficients = coplanarity(
        [b[0] - a[0], b[1] - a[1]],
        [d[0] - c[0], d[1] - c[1]],
        [c[0] - a[0], c[1] - a[1]],
    );
    earliest_contact(coefficients, travel, tolerance, |t| {
        let (first, second, _, _) =
            closest_points_on_segments(at(&a, t), at(&b, t), at(&c, t), at(&d, t));
        (first - second).norm()
    })
}

#[cfg(all(test, any(feature = "std", feature = "libm")))]
mod tests {
    use super::*;
    use azmath::Vector;

    fn v(x: f64, y: f64, z: f64) -> Vec3<f64> {
        Vector([x, y, z])
    }

    fn still(p: Vec3<f64>) -> Motion<f64> {
        [p, p]
    }

    const TOL: f64 = 1.0e-6;

    #[test]
    fn point_through_triangle() {
        let (a, b, c) = (still(v(0.0, 0.0, 0.0)), still(v(2.0, 0.0, 0.0)), still(v(0.0, 2.0, 0.0)));
        // Senkrecht durch die Mitte: Berührung nach der halben Zeit.
        let hit = vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, -1.0)], a, b, c, TOL);
        assert!(hit.is_some_and(|t| (t - 0.5).abs() < 1.0e-5));
        // Sehr schnell und weit: Die Endlage zeigt nichts, CCD findet es.
        let hit = vertex_triangle([v(0.5, 0.5, 1000.0), v(0.5, 0.5, -3000.0)], a, b, c, TOL);
        assert!(hit.is_some_and(|t| (t - 0.25).abs() < 1.0e-5));
        // Schräg, aber durch das Dreieck.
        let hit = vertex_triangle([v(0.0, 0.2, 1.0), v(1.0, 0.2, -1.0)], a, b, c, TOL);
        assert!(hit.is_some_and(|t| (t - 0.5).abs() < 1.0e-5));
        // Daneben vorbei, davor angehalten, wegbewegt.
        assert!(vertex_triangle([v(3.0, 3.0, 1.0), v(3.0, 3.0, -1.0)], a, b, c, TOL).is_none());
        assert!(vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 0.1)], a, b, c, TOL).is_none());
        assert!(vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 5.0)], a, b, c, TOL).is_none());
        // Weit weg und langsam: schneller Ausschluss.
        assert!(vertex_triangle([v(0.5, 0.5, 9.0), v(0.5, 0.5, 8.9)], a, b, c, TOL).is_none());
        // Schon am Anfang in Berührung.
        assert_eq!(vertex_triangle(still(v(0.5, 0.5, 0.0)), a, b, c, TOL), Some(0.0));
        // Genau bis zur Fläche.
        let hit = vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 0.0)], a, b, c, TOL);
        assert!(hit.is_some_and(|t| t > 0.99));
    }

    #[test]
    fn coplanar_and_common_motion() {
        let (a, b, c) = (v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0));
        let p = v(3.0, 3.0, 0.0);
        // Alles fliegt gemeinsam, schnell und in der eigenen Ebene oder quer
        // dazu: Nichts nähert sich, also keine Berührung.
        for flight in [v(500.0, -300.0, 0.0), v(0.0, 0.0, -900.0), v(70.0, 20.0, 40.0)] {
            let moved = |x: Vec3<f64>| [x, x + flight];
            assert!(vertex_triangle(moved(p), moved(a), moved(b), moved(c), TOL).is_none());
            assert!(edge_edge(moved(a), moved(b), moved(p), moved(v(4.0, 3.0, 0.0)), TOL).is_none());
        }
        // In der Ebene des Dreiecks gleitet der Punkt hinein: Das Polynom ist
        // identisch null, der Abstand entscheidet.
        let hit = vertex_triangle([p, v(0.5, 0.5, 0.0)], still(a), still(b), still(c), TOL);
        assert!(hit.is_some_and(|t| (t - 0.8).abs() < 1.0e-3));
        // Dasselbe, während alles zusammen davonfliegt.
        let flight = v(0.0, 0.0, 250.0);
        let moved = |x: Vec3<f64>| [x, x + flight];
        let hit = vertex_triangle([p, v(0.5, 0.5, 0.0) + flight], moved(a), moved(b), moved(c), TOL);
        assert!(hit.is_some_and(|t| (t - 0.8).abs() < 1.0e-3));
        // In der Ebene, aber am Dreieck vorbei.
        assert!(vertex_triangle([p, v(3.0, -3.0, 0.0)], still(a), still(b), still(c), TOL).is_none());
        // Zwei Kanten in einer Ebene laufen ineinander.
        let hit = edge_edge(still(a), still(b), [v(1.0, 1.0, 0.0), v(1.0, -1.0, 0.0)], [v(1.0, 2.0, 0.0), v(1.0, 0.0, 0.0)], TOL);
        assert!(hit.is_some_and(|t| (t - 0.5).abs() < 1.0e-3));
    }

    #[test]
    fn moving_triangle() {
        // Das Dreieck fliegt dem ruhenden Punkt entgegen.
        let p = still(v(0.5, 0.5, 0.0));
        let lift = v(0.0, 0.0, 4.0);
        let (a, b, c) = (v(0.0, 0.0, -1.0), v(2.0, 0.0, -1.0), v(0.0, 2.0, -1.0));
        let hit = vertex_triangle(p, [a, a + lift], [b, b + lift], [c, c + lift], TOL);
        assert!(hit.is_some_and(|t| (t - 0.25).abs() < 1.0e-5));
        // Dreieck dreht sich um eine Kante weg vom Punkt: keine Berührung.
        let p = still(v(0.5, 0.5, 1.0));
        let (a, b, c) = (v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0));
        let miss = vertex_triangle(p, still(a), still(b), [c, v(0.0, 2.0, -3.0)], TOL);
        assert!(miss.is_none());
        // … und zum Punkt hin: Berührung.
        let hit = vertex_triangle(p, still(a), still(b), [c, v(0.0, 2.0, 12.0)], TOL);
        assert!(hit.is_some());
    }

    #[test]
    fn edges() {
        // Eine Kante fällt quer durch eine ruhende.
        let (c, d) = (still(v(0.0, -1.0, 0.0)), still(v(0.0, 1.0, 0.0)));
        let hit = edge_edge(
            [v(-1.0, 0.0, 1.0), v(-1.0, 0.0, -1.0)],
            [v(1.0, 0.0, 1.0), v(1.0, 0.0, -1.0)],
            c,
            d,
            TOL,
        );
        assert!(hit.is_some_and(|t| (t - 0.5).abs() < 1.0e-5));
        // Sehr schnell.
        let hit = edge_edge(
            [v(-1.0, 0.0, 100.0), v(-1.0, 0.0, -300.0)],
            [v(1.0, 0.0, 100.0), v(1.0, 0.0, -300.0)],
            c,
            d,
            TOL,
        );
        assert!(hit.is_some_and(|t| (t - 0.25).abs() < 1.0e-5));
        // Seitlich versetzt: Die Geraden kreuzen sich, die Strecken nicht.
        let miss = edge_edge(
            [v(2.0, 0.0, 1.0), v(2.0, 0.0, -1.0)],
            [v(4.0, 0.0, 1.0), v(4.0, 0.0, -1.0)],
            c,
            d,
            TOL,
        );
        assert!(miss.is_none());
        // Bleibt darüber.
        let miss = edge_edge(
            [v(-1.0, 0.0, 1.0), v(-1.0, 0.0, 0.5)],
            [v(1.0, 0.0, 1.0), v(1.0, 0.0, 0.5)],
            c,
            d,
            TOL,
        );
        assert!(miss.is_none());
        // Schon in Berührung.
        assert_eq!(
            edge_edge(still(v(-1.0, 0.0, 0.0)), still(v(1.0, 0.0, 0.0)), c, d, TOL),
            Some(0.0)
        );
    }

    /// Viele zufällige Bewegungen: Wo die Endlage das Dreieck durchstoßen
    /// hat, muss CCD eine Berührung gemeldet haben.
    #[test]
    fn never_misses_a_crossing() {
        let (a, b, c) = (v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0));
        let mut seed = 12345_u64;
        let mut random = || {
            // Einfacher Zufallsgenerator (LCG), Werte in [−1, 1).
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            ((seed >> 11) as f64 / (1_u64 << 53) as f64) * 2.0 - 1.0
        };
        let mut crossings = 0;
        for _ in 0..5000 {
            let start = v(random(), random(), random() * 0.5 + 1.0);
            let end = v(random(), random(), random() * 0.5 - 1.0);
            let crosses = crate::geometry::segment_intersects_triangle(start, end, a, b, c);
            let detected = vertex_triangle([start, end], still(a), still(b), still(c), TOL).is_some();
            if crosses {
                crossings += 1;
                assert!(detected);
            }
        }
        assert!(crossings > 100);
    }
}
