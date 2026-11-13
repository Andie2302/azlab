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
//! Bei positivem Kontaktabstand müssen die Teile nicht in einer Ebene
//! liegen. Deshalb wird der Abstand während des gesamten Schritts verfolgt:
//! Es wird jeweils höchstens so weit vorgerückt, wie es die obere Schranke
//! der relativen Bewegung erlaubt. Das gilt auch für koplanare Bewegungen.
//!
//! Gerechnet wird mit der Bewegung der Teile gegeneinander. Eine gemeinsame
//! Verschiebung – das Tuch fliegt als Ganzes – zählt nicht als Annäherung.
//!
//! Das Ergebnis ([`Contact`]) kennt drei Fälle: sicher frei, sicher berührt
//! oder nicht entschieden, wenn die Höchstzahl der Schritte erreicht ist.
//! Die Prüfung meldet nie „frei“, wenn sich die Teile durchdringen.

use crate::geometry::{Vec3, closest_points_on_segments, point_triangle_distance, segment_segment_distance};
use azmath::scalar::RealElementary;

/// Lage am Anfang und am Ende des Zeitschritts.
pub type Motion<T> = [Vec3<T>; 2];

/// Ergebnis der kontinuierlichen Prüfung eines Paars.
///
/// Drei Fälle, weil die Prüfung nach einer festen Zahl von Schritten endet:
/// Sie kann eine Berührung sicher ausschließen, sicher feststellen – oder
/// keins von beiden, wenn sich die Teile schleichend nähern.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Contact<T> {
    /// Sicher keine Berührung im ganzen Schritt.
    Free,
    /// Sicher eine Berührung: Zu diesem Zeitpunkt in [0, 1] ist der Abstand
    /// erstmals höchstens die Toleranz.
    Touch(T),
    /// Nicht entschieden: Bis zu diesem Zeitpunkt in [0, 1] berühren sich
    /// die Teile sicher nicht; für die Zeit danach hat die Prüfung ihre
    /// Höchstzahl an Schritten erreicht, ohne eine Berührung zu finden oder
    /// auszuschließen. Auch ein nicht berechenbarer Abstand (NaN) endet so.
    Undecided(T),
}

impl<T: Copy> Contact<T> {
    /// `true`, wenn eine Berührung sicher ausgeschlossen ist.
    #[inline]
    pub fn is_free(&self) -> bool {
        matches!(self, Self::Free)
    }

    /// Zeitpunkt der sicher festgestellten Berührung.
    #[inline]
    pub fn touch(&self) -> Option<T> {
        match *self {
            Self::Touch(time) => Some(time),
            _ => None,
        }
    }

    /// Zeitpunkt, bis zu dem die Bewegung sicher frei ist: der Zeitpunkt
    /// der Berührung oder der, an dem die Prüfung abbrach. `None`, wenn der
    /// ganze Schritt frei ist.
    ///
    /// Wer nichts durchdringen lassen will, behandelt beide Fälle gleich
    /// und beendet die Bewegung vor diesem Zeitpunkt.
    #[inline]
    pub fn safe_until(&self) -> Option<T> {
        match *self {
            Self::Free => None,
            Self::Touch(time) | Self::Undecided(time) => Some(time),
        }
    }
}

/// Lage zur Zeit `t` in [0, 1].
#[inline]
fn at<T: RealElementary>(motion: &Motion<T>, t: T) -> Vec3<T> {
    motion[0] + (motion[1] - motion[0]) * t
}

/// Größte Strecke, die einer der Punkte im Schritt zurücklegt.
fn reach<T: RealElementary>(motions: &[&Motion<T>]) -> T {
    motions.iter().fold(T::ZERO, |longest, motion| longest.max((motion[1] - motion[0]).norm()))
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

/// Höchstzahl der Schritte beim Vorrücken.
const MAX_ADVANCES: usize = 64;

/// Sucht in [`from`, `to`] den frühesten Zeitpunkt, an dem `distance(t)`
/// höchstens `tolerance` beträgt, durch vorsichtiges Vorrücken.
///
/// Der Abstand ändert sich je Zeiteinheit höchstens um `reach`. Vom
/// Zeitpunkt `t` aus kann die Berührung also frühestens nach
/// `(Abstand − Toleranz) / reach` eintreten; so weit wird vorgerückt. Das
/// übersieht keine Berührung, für keine Art von Bewegung.
///
/// Als Berührung gilt auch, wenn der Abstand die Toleranz bis auf Rundung
/// erreicht hat und ein weiterer Schritt die Zeit nicht mehr merklich
/// voranbrächte. Nähern sich die Teile dagegen schleichend, werden die
/// Schritte nur allmählich kleiner; nach [`MAX_ADVANCES`] Schritten endet
/// die Suche mit [`Contact::Undecided`].
fn advance<T: RealElementary>(
    from: T,
    to: T,
    reach: T,
    tolerance: T,
    distance: impl Fn(T) -> T,
) -> Contact<T> {
    let mut t = from;
    for _ in 0..MAX_ADVANCES {
        let gap = distance(t) - tolerance;
        if gap.is_nan() {
            return Contact::Undecided(t);
        }
        if gap <= T::ZERO {
            return Contact::Touch(t);
        }
        if reach <= T::ZERO {
            return Contact::Free;
        }
        let step = gap / reach;
        if t + step > to {
            return Contact::Free;
        }
        // Der Abstand liegt nur noch um Rundung über der Toleranz: Weiter
        // lässt sich die Zeit nicht auflösen, die Berührung ist erreicht.
        // (Ohne diese Prüfung endete auch ein gerader Aufprall unentschieden,
        // weil der Abstand die Toleranz rechnerisch nie genau trifft.)
        if step <= T::EPSILON * T::from_f64(16.0) {
            return Contact::Touch(t);
        }
        t += step;
    }
    Contact::Undecided(t)
}

/// Wann berührt der bewegte Punkt `p` das bewegte Dreieck `a`, `b`, `c`?
///
/// Jedes Argument ist die Lage am Anfang und am Ende des Schritts.
/// `tolerance` ist der Abstand, der als Berührung gilt (größer als null,
/// klein gegen die Abmessungen). Das Ergebnis unterscheidet drei Fälle,
/// siehe [`Contact`]: sicher frei, sicher berührt (mit dem frühesten
/// Zeitpunkt) oder nicht entschieden.
pub fn vertex_triangle<T: RealElementary>(
    p: Motion<T>,
    a: Motion<T>,
    b: Motion<T>,
    c: Motion<T>,
    tolerance: T,
) -> Contact<T> {
    let [p, a, b, c] = relative([p, a, b, c]);
    // Schneller Ausschluss: Der Abstand kann höchstens um die Wege schrumpfen,
    // die Punkt und Dreieck zurücklegen.
    let travel = (p[1] - p[0]).norm() + reach(&[&a, &b, &c]);
    let start = point_triangle_distance(p[0], a[0], b[0], c[0]);
    if start > travel + tolerance {
        return Contact::Free;
    }
    if start <= tolerance {
        return Contact::Touch(T::ZERO);
    }
    advance(T::ZERO, T::ONE, travel, tolerance, |t| {
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
) -> Contact<T> {
    let [a, b, c, d] = relative([a, b, c, d]);
    let travel = reach(&[&a, &b]) + reach(&[&c, &d]);
    let start = segment_segment_distance(a[0], b[0], c[0], d[0]);
    if start > travel + tolerance {
        return Contact::Free;
    }
    if start <= tolerance {
        return Contact::Touch(T::ZERO);
    }
    advance(T::ZERO, T::ONE, travel, tolerance, |t| {
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
        assert!(hit.touch().is_some_and(|t| (t - 0.5).abs() < 1.0e-5));
        // Sehr schnell und weit: Die Endlage zeigt nichts, CCD findet es.
        let hit = vertex_triangle([v(0.5, 0.5, 1000.0), v(0.5, 0.5, -3000.0)], a, b, c, TOL);
        assert!(hit.touch().is_some_and(|t| (t - 0.25).abs() < 1.0e-5));
        // Schräg, aber durch das Dreieck.
        let hit = vertex_triangle([v(0.0, 0.2, 1.0), v(1.0, 0.2, -1.0)], a, b, c, TOL);
        assert!(hit.touch().is_some_and(|t| (t - 0.5).abs() < 1.0e-5));
        // Daneben vorbei, davor angehalten, wegbewegt.
        assert!(vertex_triangle([v(3.0, 3.0, 1.0), v(3.0, 3.0, -1.0)], a, b, c, TOL).is_free());
        assert!(vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 0.1)], a, b, c, TOL).is_free());
        assert!(vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 5.0)], a, b, c, TOL).is_free());
        // Weit weg und langsam: schneller Ausschluss.
        assert!(vertex_triangle([v(0.5, 0.5, 9.0), v(0.5, 0.5, 8.9)], a, b, c, TOL).is_free());
        // Schon am Anfang in Berührung.
        assert_eq!(vertex_triangle(still(v(0.5, 0.5, 0.0)), a, b, c, TOL), Contact::Touch(0.0));
        // Genau bis zur Fläche.
        let hit = vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 0.0)], a, b, c, TOL);
        assert!(hit.touch().is_some_and(|t| t > 0.99));
    }

    #[test]
    fn coplanar_and_common_motion() {
        let (a, b, c) = (v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0));
        let p = v(3.0, 3.0, 0.0);
        // Alles fliegt gemeinsam, schnell und in der eigenen Ebene oder quer
        // dazu: Nichts nähert sich, also keine Berührung.
        for flight in [v(500.0, -300.0, 0.0), v(0.0, 0.0, -900.0), v(70.0, 20.0, 40.0)] {
            let moved = |x: Vec3<f64>| [x, x + flight];
            assert!(vertex_triangle(moved(p), moved(a), moved(b), moved(c), TOL).is_free());
            assert!(edge_edge(moved(a), moved(b), moved(p), moved(v(4.0, 3.0, 0.0)), TOL).is_free());
        }
        // In der Ebene des Dreiecks gleitet der Punkt hinein: Das Polynom ist
        // identisch null, der Abstand entscheidet.
        let hit = vertex_triangle([p, v(0.5, 0.5, 0.0)], still(a), still(b), still(c), TOL);
        assert!(hit.touch().is_some_and(|t| (t - 0.8).abs() < 1.0e-3));
        // Dasselbe, während alles zusammen davonfliegt.
        let flight = v(0.0, 0.0, 250.0);
        let moved = |x: Vec3<f64>| [x, x + flight];
        let hit = vertex_triangle([p, v(0.5, 0.5, 0.0) + flight], moved(a), moved(b), moved(c), TOL);
        assert!(hit.touch().is_some_and(|t| (t - 0.8).abs() < 1.0e-3));
        // In der Ebene, aber am Dreieck vorbei.
        assert!(vertex_triangle([p, v(3.0, -3.0, 0.0)], still(a), still(b), still(c), TOL).is_free());
        // Zwei Kanten in einer Ebene laufen ineinander.
        let hit = edge_edge(still(a), still(b), [v(1.0, 1.0, 0.0), v(1.0, -1.0, 0.0)], [v(1.0, 2.0, 0.0), v(1.0, 0.0, 0.0)], TOL);
        assert!(hit.touch().is_some_and(|t| (t - 0.5).abs() < 1.0e-3));
    }

    #[test]
    fn moving_triangle() {
        // Das Dreieck fliegt dem ruhenden Punkt entgegen.
        let p = still(v(0.5, 0.5, 0.0));
        let lift = v(0.0, 0.0, 4.0);
        let (a, b, c) = (v(0.0, 0.0, -1.0), v(2.0, 0.0, -1.0), v(0.0, 2.0, -1.0));
        let hit = vertex_triangle(p, [a, a + lift], [b, b + lift], [c, c + lift], TOL);
        assert!(hit.touch().is_some_and(|t| (t - 0.25).abs() < 1.0e-5));
        // Dreieck dreht sich um eine Kante weg vom Punkt: keine Berührung.
        let p = still(v(0.5, 0.5, 1.0));
        let (a, b, c) = (v(0.0, 0.0, 0.0), v(2.0, 0.0, 0.0), v(0.0, 2.0, 0.0));
        let miss = vertex_triangle(p, still(a), still(b), [c, v(0.0, 2.0, -3.0)], TOL);
        assert!(miss.is_free());
        // … und zum Punkt hin: Berührung.
        let hit = vertex_triangle(p, still(a), still(b), [c, v(0.0, 2.0, 12.0)], TOL);
        assert!(!hit.is_free());
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
        assert!(hit.touch().is_some_and(|t| (t - 0.5).abs() < 1.0e-5));
        // Sehr schnell.
        let hit = edge_edge(
            [v(-1.0, 0.0, 100.0), v(-1.0, 0.0, -300.0)],
            [v(1.0, 0.0, 100.0), v(1.0, 0.0, -300.0)],
            c,
            d,
            TOL,
        );
        assert!(hit.touch().is_some_and(|t| (t - 0.25).abs() < 1.0e-5));
        // Seitlich versetzt: Die Geraden kreuzen sich, die Strecken nicht.
        let miss = edge_edge(
            [v(2.0, 0.0, 1.0), v(2.0, 0.0, -1.0)],
            [v(4.0, 0.0, 1.0), v(4.0, 0.0, -1.0)],
            c,
            d,
            TOL,
        );
        assert!(miss.is_free());
        // Bleibt darüber.
        let miss = edge_edge(
            [v(-1.0, 0.0, 1.0), v(-1.0, 0.0, 0.5)],
            [v(1.0, 0.0, 1.0), v(1.0, 0.0, 0.5)],
            c,
            d,
            TOL,
        );
        assert!(miss.is_free());
        // Schon in Berührung.
        assert_eq!(
            edge_edge(still(v(-1.0, 0.0, 0.0)), still(v(1.0, 0.0, 0.0)), c, d, TOL),
            Contact::Touch(0.0)
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
            let detected = !vertex_triangle([start, end], still(a), still(b), still(c), TOL).is_free();
            if crosses {
                crossings += 1;
                assert!(detected);
            }
        }
        assert!(crossings > 100);
    }
    #[test]
    fn positive_tolerance_without_plane_crossing() {
        let (a, b, c) = (still(v(0.0, 0.0, 0.0)), still(v(2.0, 0.0, 0.0)), still(v(0.0, 2.0, 0.0)));
        let tolerance = 0.01;
        // Erreicht den Kontaktabstand, ohne die Dreiecksebene zu erreichen.
        let hit = vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 0.005)], a, b, c, tolerance);
        let expected = (1.0 - tolerance) / (1.0 - 0.005);
        assert!(hit.touch().is_some_and(|t| (t - expected).abs() < 1.0e-6));

        // Kontakt nur mitten im Schritt: Die Endlagen allein reichen nicht.
        let hit = vertex_triangle([v(-1.0, 0.5, 0.005), v(3.0, 0.5, 0.005)], a, b, c, tolerance);
        assert!(hit.touch().is_some_and(|t| t > 0.24 && t < 0.26));
        // Dieselbe Bahn mit ausreichendem Abstand bleibt frei.
        assert!(vertex_triangle([v(-1.0, 0.5, 0.2), v(3.0, 0.5, 0.2)], a, b, c, tolerance).is_free());
    }

    #[test]
    fn edge_tolerance_without_plane_crossing() {
        let (c, d) = (still(v(0.0, -1.0, 0.0)), still(v(0.0, 1.0, 0.0)));
        let tolerance = 0.01;
        let hit = edge_edge(
            [v(-1.0, 0.0, 1.0), v(-1.0, 0.0, 0.005)],
            [v(1.0, 0.0, 1.0), v(1.0, 0.0, 0.005)],
            c, d, tolerance,
        );
        let expected = (1.0 - tolerance) / (1.0 - 0.005);
        assert!(hit.touch().is_some_and(|t| (t - expected).abs() < 1.0e-6));
        // Windschiefe Kanten erreichen den Kontaktabstand nur mitten im Schritt.
        let hit = edge_edge(
            [v(-1.0, -2.0, 0.005), v(-1.0, 2.0, 0.005)],
            [v(1.0, -2.0, 0.005), v(1.0, 2.0, 0.005)],
            c, d, tolerance,
        );
        assert!(hit.touch().is_some_and(|t| t > 0.24 && t < 0.26));
        assert!(edge_edge(
            [v(-1.0, -2.0, 0.2), v(-1.0, 2.0, 0.2)],
            [v(1.0, -2.0, 0.2), v(1.0, 2.0, 0.2)],
            c, d, tolerance,
        ).is_free());
    }

    #[test]
    fn three_outcomes() {
        let (a, b, c) = (still(v(0.0, 0.0, 0.0)), still(v(2.0, 0.0, 0.0)), still(v(0.0, 2.0, 0.0)));
        let tolerance = 0.01;
        // Sicher frei, sicher berührt.
        let free = vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, 0.5)], a, b, c, tolerance);
        assert_eq!(free, Contact::Free);
        assert!(free.is_free() && free.touch().is_none() && free.safe_until().is_none());
        let touch = vertex_triangle([v(0.5, 0.5, 1.0), v(0.5, 0.5, -1.0)], a, b, c, tolerance);
        assert!(matches!(touch, Contact::Touch(t) if (t - 0.495).abs() < 1.0e-6));
        assert_eq!(touch.touch(), touch.safe_until());

        // Schleichend: Der Punkt gleitet weit und dicht über das Dreieck,
        // knapp außerhalb der Toleranz. Jeder Schritt darf nur winzig
        // vorrücken; die Prüfung endet unentschieden.
        let height = tolerance * (1.0 + 1.0e-9);
        let glide = [v(0.1, 0.1, height), v(0.9, 0.9, height)];
        let undecided = vertex_triangle(glide, a, b, c, tolerance);
        assert!(matches!(undecided, Contact::Undecided(t) if (0.0..1.0).contains(&t)));
        assert!(!undecided.is_free() && undecided.touch().is_none());
        // Bis zum gemeldeten Zeitpunkt ist die Bahn nachweislich frei.
        assert!(undecided.safe_until().is_some_and(|time| {
            (0..=10).all(|step| {
                let t = time * f64::from(step) / 10.0;
                point_triangle_distance(at(&glide, t), a[0], b[0], c[0]) > tolerance
            })
        }));
        // Dasselbe für zwei Kanten.
        let (p, q) = (still(v(0.0, -1.0, 0.0)), still(v(0.0, 1.0, 0.0)));
        let undecided = edge_edge(
            [v(-1.0, -0.5, height), v(-1.0, 0.5, height)],
            [v(1.0, -0.5, height), v(1.0, 0.5, height)],
            p, q, tolerance,
        );
        assert!(matches!(undecided, Contact::Undecided(_)));
        // Nicht berechenbarer Abstand: ebenfalls unentschieden, nie „frei“.
        let broken = vertex_triangle([v(0.5, 0.5, 1.0), v(f64::NAN, 0.5, -1.0)], a, b, c, tolerance);
        assert!(!broken.is_free());
    }
}
