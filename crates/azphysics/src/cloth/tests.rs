//! Tests des Gewebes: Ruhelage, Hindernisse und vor allem die Zusage, dass
//! nichts durchdrungen wird.

use super::{Cloth, ClothParams, Collider};
use azmath::Vector;
use azmath::scalar::RealElementary;

fn v<T: RealElementary>(x: f64, y: f64, z: f64) -> Vector<T, 3> {
    Vector([T::from_f64(x), T::from_f64(y), T::from_f64(z)])
}

/// Tuch in der Ebene y = `height`, 1 m × 1 m, `n × n` Teilchen.
fn sheet<T: RealElementary>(n: usize, height: f64) -> Cloth<T> {
    let spacing = 1.0 / (n - 1) as f64;
    Cloth::grid(n, n, v(-0.5, height, -0.5), v(spacing, 0.0, 0.0), v(0.0, 0.0, spacing), T::from_f64(0.01))
}

fn hanging<T: RealElementary>() {
    let f = T::from_f64;
    // An der oberen Kante aufgehängt, hängt senkrecht nach unten.
    let n = 9;
    let spacing = 1.0 / (n - 1) as f64;
    let mut cloth: Cloth<T> =
        Cloth::grid(n, n, v(-0.5, 1.0, 0.0), v(spacing, 0.0, 0.0), v(0.0, 0.0, spacing), f(0.01));
    for column in 0..n {
        let index = cloth.index(column, 0);
        cloth.pin(index);
    }
    assert_eq!(cloth.positions().len(), n * n);
    assert_eq!(cloth.cloth_triangle_count(), 2 * (n - 1) * (n - 1));
    assert_eq!(cloth.triangles().len(), cloth.cloth_triangle_count());

    for _ in 0..180 {
        cloth.step(f(1.0 / 60.0), &[]);
    }
    // Aufgehängte Teilchen bleiben, der Rest ist nach unten geschwungen.
    assert_eq!(cloth.positions()[cloth.index(0, 0)], v(-0.5, 1.0, 0.0));
    let lowest = cloth.positions().iter().fold(T::INFINITY, |low, p| low.min(p[1]));
    assert!(lowest < f(0.3));
    // Kaum gedehnt, nichts entgleist, nichts durchdrungen.
    assert!(cloth.max_strain() < f(0.1));
    assert!(cloth.positions().iter().all(|p| p.iter().all(|c| c.is_finite())));
    assert!(cloth.velocities().iter().all(|u| u.norm() < f(5.0)));
    assert_eq!(cloth.intersections(), 0);
}

fn drape<T: RealElementary>() {
    let f = T::from_f64;
    // Tuch fällt auf eine Kugel über dem Boden.
    let mut cloth = sheet::<T>(11, 1.0);
    let (center, radius) = (v::<T>(0.0, 0.5, 0.0), f(0.3));
    let colliders = [
        Collider::Sphere { center, radius },
        Collider::Plane { point: v(0.0, 0.0, 0.0), normal: v(0.0, 1.0, 0.0) },
    ];
    let mut repulsions = 0;
    for _ in 0..150 {
        repulsions += cloth.step(f(1.0 / 60.0), &colliders).repulsions;
    }
    let slack = cloth.params.thickness * f(0.5);
    for p in cloth.positions() {
        assert!((*p - center).norm() >= radius - slack);
        assert!(p[1] >= -slack);
    }
    // Die Mitte liegt auf der Kugel, die Ränder hängen tiefer.
    let top = cloth.positions()[cloth.index(5, 5)];
    assert!(top[1] > f(0.75) && top[1] < f(0.9));
    assert!(cloth.positions()[cloth.index(0, 0)][1] < top[1]);
    assert!(cloth.max_strain() < f(0.2));
    assert_eq!(cloth.intersections(), 0);
    let _ = repulsions;
}

/// Kern der Sache: Ein Tuch wird mit sehr hoher Geschwindigkeit auf ein
/// dünnes Hindernis geschossen. Je Teilschritt legt es ein Vielfaches seiner
/// Dicke zurück.
fn tunnelling<T: RealElementary>() {
    let f = T::from_f64;
    let launch = |ccd: bool| {
        let mut cloth = sheet::<T>(5, 0.5);
        cloth.params = ClothParams { ccd, gravity: v(0.0, 0.0, 0.0), damping: T::ZERO, ..ClothParams::default() };
        // Großes Dreieck in der Ebene y = 0, weit über das Tuch hinaus.
        cloth.add_obstacle_triangle(v(-5.0, 0.0, -5.0), v(5.0, 0.0, -5.0), v(0.0, 0.0, 8.0));
        cloth.set_all_velocities(v(0.0, -300.0, 0.0));
        let mut stats = (0, 0);
        for _ in 0..30 {
            let step = cloth.step(f(1.0 / 60.0), &[]);
            stats = (stats.0 + step.ccd_hits, stats.1 + step.full_rollbacks);
        }
        (cloth, stats)
    };

    // Ohne Sicherung: Das Tuch ist zwischen zwei Teilschritten einfach
    // hindurchgesprungen.
    let (unsafe_cloth, _) = launch(false);
    let moving = 25;
    assert!(unsafe_cloth.positions()[..moving].iter().all(|p| p[1] < T::ZERO));

    // Mit Sicherung: Kein Teilchen kommt je auf die andere Seite. Das Tuch
    // ist dabei wirklich bis an das Hindernis herangeflogen und nicht etwa
    // in der Luft stehen geblieben.
    let (safe_cloth, (hits, _)) = launch(true);
    let (low, high) = safe_cloth.positions()[..moving]
        .iter()
        .fold((T::INFINITY, T::NEG_INFINITY), |(low, high), p| (low.min(p[1]), high.max(p[1])));
    assert!(low > T::ZERO && high < f(0.1), "Höhen {:?} bis {:?}", low.to_f64(), high.to_f64());
    assert!(hits > 0);

    // Ein Tuch im freien Flug stößt nicht mit sich selbst zusammen: Die
    // Sicherung greift nicht ein, es kommt genau so weit wie ohne sie.
    let fly = |ccd: bool| {
        let mut cloth = sheet::<T>(5, 0.5);
        cloth.params = ClothParams { ccd, gravity: v(0.0, 0.0, 0.0), damping: T::ZERO, ..ClothParams::default() };
        cloth.set_all_velocities(v(30.0, -300.0, 10.0));
        let hits = (0..6).fold(0, |hits, _| hits + cloth.step(f(1.0 / 60.0), &[]).ccd_hits);
        (cloth.positions()[12], hits)
    };
    let ((free, _), (checked, hits)) = (fly(false), fly(true));
    assert_eq!(hits, 0);
    assert!((free - checked).norm() < f(1.0e-3) && free[1] < f(-20.0));
    assert_eq!(safe_cloth.intersections(), 0);
    assert_eq!(safe_cloth.triangles().len(), safe_cloth.cloth_triangle_count() + 1);
}

/// Das Tuch wird in sich selbst gefaltet: eine Hälfte fliegt schnell durch
/// die Ebene der anderen.
fn folding<T: RealElementary>() {
    let f = T::from_f64;
    let n = 9;
    let mut cloth = sheet::<T>(n, 0.5);
    cloth.params.gravity = v(0.0, 0.0, 0.0);
    // Linke Hälfte festhalten, rechte Hälfte um die Mittellinie nach oben
    // und zurück über die linke klappen.
    for row in 0..n {
        for column in 0..n {
            let index = cloth.index(column, row);
            if column < n / 2 {
                cloth.pin(index);
            } else {
                cloth.set_velocity(index, v(-40.0, 25.0, 0.0));
            }
        }
    }
    let mut hits = 0;
    for _ in 0..60 {
        hits += cloth.step(f(1.0 / 60.0), &[]).ccd_hits;
        assert_eq!(cloth.intersections(), 0);
    }
    assert!(cloth.positions().iter().all(|p| p.iter().all(|c| c.is_finite())));
    let _ = hits;
}

fn editing<T: RealElementary>() {
    let f = T::from_f64;
    let mut cloth = sheet::<T>(3, 1.0);
    // Zu kleine Gitter werden auf 2 × 2 angehoben.
    let tiny: Cloth<T> = Cloth::grid(1, 1, v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 0.0, 1.0), f(1.0));
    assert_eq!(tiny.positions().len(), 4);
    // Festgehaltene Teilchen nehmen keine Geschwindigkeit an, lassen sich
    // aber versetzen.
    cloth.pin(0);
    cloth.set_velocity(0, v(1.0, 2.0, 3.0));
    cloth.set_velocity(1, v(1.0, 2.0, 3.0));
    cloth.set_velocity(99, v(1.0, 2.0, 3.0));
    assert_eq!(cloth.velocities()[0], v(0.0, 0.0, 0.0));
    assert_eq!(cloth.velocities()[1], v(1.0, 2.0, 3.0));
    cloth.set_position(0, v(9.0, 9.0, 9.0));
    cloth.set_position(99, v(9.0, 9.0, 9.0));
    cloth.pin(99);
    assert_eq!(cloth.positions()[0], v(9.0, 9.0, 9.0));
    // Ungültige Schrittweite ändert nichts.
    let before = cloth.positions()[1];
    cloth.step(T::ZERO, &[]);
    cloth.step(T::NAN, &[]);
    assert_eq!(cloth.positions()[1], before);
    assert_eq!(cloth.max_strain(), cloth.max_strain());
}

#[test]
fn cloth_f32() {
    hanging::<f32>();
    drape::<f32>();
    tunnelling::<f32>();
    folding::<f32>();
    editing::<f32>();
}

#[test]
fn cloth_f64() {
    hanging::<f64>();
    drape::<f64>();
    tunnelling::<f64>();
    folding::<f64>();
    editing::<f64>();
}

#[test]
fn observing_does_not_change_cloth() {
    let mut plain = sheet::<f64>(4, 0.03);
    let mut observed = plain.clone();
    let colliders = [Collider::Plane { point: v(0.0, 0.0, 0.0), normal: v(0.0, 1.0, 0.0) }];
    let mut phases = alloc::vec::Vec::new();
    for _ in 0..10 {
        let a = plain.step(1.0 / 60.0, &colliders);
        let b = observed.step_observed(1.0 / 60.0, &colliders, &mut |phase| phases.push(phase));
        assert_eq!(a, b);
        assert_eq!(plain.positions(), observed.positions());
        assert_eq!(plain.velocities(), observed.velocities());
    }
    assert!(phases.contains(&super::StepPhase::Proximity));
    assert!(phases.contains(&super::StepPhase::Ccd));
    phases.clear();
    observed.step_observed(0.0, &colliders, &mut |phase| phases.push(phase));
    assert!(phases.is_empty());
}

#[test]
fn motion_restore_validates_before_modifying() {
    let mut cloth = sheet::<f64>(4, 1.0);
    let positions = cloth.positions().to_vec();
    let velocities = cloth.velocities().to_vec();
    assert_eq!(cloth.restore_motion(&positions[..2], &velocities), Err(crate::MotionError::DimensionMismatch));
    let mut bad = positions.clone();
    bad[0][0] = f64::NAN;
    assert_eq!(cloth.restore_motion(&bad, &velocities), Err(crate::MotionError::NonFinite));
    assert_eq!(cloth.positions(), positions);
    cloth.step(1.0 / 60.0, &[]);
    assert!(cloth.restore_motion(&positions, &velocities).is_ok());
    assert_eq!(cloth.positions(), positions);
    assert_eq!(cloth.velocities(), velocities);
    let mut expected = sheet::<f64>(4, 1.0);
    for _ in 0..3 {
        cloth.step(1.0 / 60.0, &[]);
        expected.step(1.0 / 60.0, &[]);
    }
    assert_eq!(cloth.positions(), expected.positions());
}


/// Der Diagnosemodus muss echte freie Teilchen liefern: Ein gedehntes Netz
/// darf bei null Iterationen keine versteckte Mindestiteration ausführen.
/// Jedes freie Teilchen behält seine Geschwindigkeit bitgenau und liegt bis
/// auf Rundung bei Anfangslage plus Geschwindigkeit mal Zeit.
fn zero_iterations_preserves_free_particle_path<T: RealElementary>() {
    let mut cloth = sheet::<T>(3, 0.0);
    cloth.params.iterations = 0;
    cloth.params.ccd = false;
    cloth.params.self_collision = false;
    cloth.params.gravity = v(0.0, 0.0, 0.0);
    cloth.params.damping = T::ZERO;
    cloth.pin(0);
    cloth.set_position(4, v(0.7, 0.3, -0.2));
    let velocity = v::<T>(0.15, 0.0, 0.05);
    cloth.set_all_velocities(velocity);
    let start = cloth.positions().to_vec();
    let dt = T::from_f64(1.0 / 60.0);
    for step in 1..=3 {
        cloth.step(dt, &[]);
        let time = dt * T::from_f64(f64::from(step));
        assert_eq!(cloth.positions()[0], start[0]);
        assert_eq!(cloth.velocities()[0], v(0.0, 0.0, 0.0));
        for (index, &began) in start.iter().enumerate().skip(1) {
            assert_eq!(cloth.velocities()[index], velocity);
            let error = cloth.positions()[index] - (began + velocity * time);
            assert!(error.norm() <= T::EPSILON * T::from_f64(4.0), "Teilchen {index}: {:?}", error.norm().to_f64());
        }
    }
}

/// Freier Flug weit weg vom Ursprung, mit allen Korrekturen eingeschaltet.
///
/// Dort ist das Raster der Lage grob. Früher wurde die Geschwindigkeit in
/// jedem Teilschritt aus dem gerundeten Weg zurückgerechnet und rastete auf
/// dieses Raster ein; Nachbarn wurden verschieden gerundet, und die
/// Verbindungen verteilten Stöße in Rastergröße. Das Tuch trieb stetig ab –
/// in jeder Richtung, in der es sich bewegte, und umso stärker, je weiter
/// es vom Ursprung entfernt war.
fn free_flight_does_not_drift<T: RealElementary>() {
    let f = T::from_f64;
    let n = 8;
    let spacing = 1.0 / (n - 1) as f64;
    let origin = v::<T>(250.0, -40.0, 90.0);
    let mut cloth: Cloth<T> = Cloth::grid(n, n, origin, v(spacing, 0.0, 0.0), v(0.0, 0.0, spacing), f(0.005));
    cloth.params.gravity = v(0.0, 0.0, 0.0);
    cloth.params.damping = T::ZERO;
    let velocity = v::<T>(0.15, -0.02, 0.05);
    cloth.set_all_velocities(velocity);
    let start = cloth.positions().to_vec();
    // Ein Rasterschritt der Lage an dieser Stelle.
    let grid = f(256.0) * T::EPSILON;
    let dt = f(1.0 / 60.0);
    let steps = 240;
    let mut hits = 0;
    for _ in 0..steps {
        hits += cloth.step(dt, &[]).ccd_hits;
    }
    let time = dt * f(f64::from(steps));
    assert_eq!(hits, 0);
    for (index, (&position, &began)) in cloth.positions().iter().zip(&start).enumerate() {
        // Die Geschwindigkeit ist unverändert, bis auf das letzte Bit.
        assert_eq!(cloth.velocities()[index], velocity, "Teilchen {index}");
        // Die Lage stimmt auf wenige Rasterschritte. Die alte Rechnung lag
        // hier nach vier Sekunden um mehrere hundert daneben.
        let error = position - (began + velocity * time);
        assert!(error.norm() <= grid * f(4.0), "Teilchen {index}: {:?} Rasterschritte", (error.norm() / grid).to_f64());
    }
    assert!(cloth.max_strain() <= grid * f(64.0));
}

/// Ein feines Tuch hängt senkrecht an seiner oberen Kante, genau in seiner
/// Ruhelage. Es muss ruhig hängen bleiben.
///
/// Solange die Verbindungen stets in derselben Reihenfolge durchlaufen
/// wurden, schaukelte es sich aus der Ruhe auf: nach einer halben Sekunde
/// mehrere Meter je Sekunde, dann faltete es sich in der eigenen Ebene
/// zusammen (Dehnungen über 60 %) oder lief bis zum Überlauf davon. Mit 12
/// Teilchen je Kante trat das nicht auf, ab 24 mit den Standardwerten immer.
fn hanging_cloth_stays_at_rest<T: RealElementary>(iterations: u32) {
    let f = T::from_f64;
    let n = 24;
    let spacing = 1.0 / (n - 1) as f64;
    let mut cloth: Cloth<T> =
        Cloth::grid(n, n, v(-0.5, 1.5, 0.0), v(spacing, 0.0, 0.0), v(0.0, -spacing, 0.0), f(0.005));
    cloth.params.iterations = iterations;
    // Für diese Frage ohne Belang und in der Prüfung nur langsam.
    cloth.params.self_collision = false;
    cloth.params.ccd = false;
    for column in 0..n {
        let index = cloth.index(column, 0);
        cloth.pin(index);
    }
    let mut fastest = T::ZERO;
    // Die alte Rechnung war nach 40 Bildern bei mehreren Metern je Sekunde.
    for _ in 0..60 {
        cloth.step(f(1.0 / 60.0), &[]);
        fastest = cloth.velocities().iter().fold(fastest, |fastest, velocity| fastest.max(velocity.norm()));
        assert!(fastest.is_finite());
    }
    // Es sackt beim Loslassen etwas durch und schwingt aus; schnell wird es
    // dabei nie.
    assert!(fastest < f(0.5), "bis {:?} m/s schnell", fastest.to_f64());
    let speed = cloth.velocities().iter().fold(T::ZERO, |fastest, velocity| fastest.max(velocity.norm()));
    assert!(speed < f(0.1), "am Ende noch {:?} m/s", speed.to_f64());
    // Mit so wenigen Durchläufen gibt das Tuch etwas nach, aber nur wenig.
    assert!(cloth.max_strain() < f(0.03), "Dehnung {:?}", cloth.max_strain().to_f64());
    let lowest = cloth.positions().iter().fold(T::INFINITY, |low, point| low.min(point[1]));
    assert!(lowest > f(0.49) && lowest < f(0.5));
    // Es ist in seiner Ebene geblieben.
    assert!(cloth.positions().iter().all(|point| point[2] == T::ZERO));
}

#[test]
fn hanging_rest_f32() {
    // Eins wird auf ein Paar aufgerundet; ein einzelner Durchlauf wäre
    // nicht stabil.
    hanging_cloth_stays_at_rest::<f32>(1);
}

#[test]
fn hanging_rest_f64() {
    hanging_cloth_stays_at_rest::<f64>(2);
}

#[test]
fn free_flight_f32() {
    free_flight_does_not_drift::<f32>();
}

#[test]
fn free_flight_f64() {
    free_flight_does_not_drift::<f64>();
}

#[test]
fn zero_iterations_free_path_f32() {
    zero_iterations_preserves_free_particle_path::<f32>();
}

#[test]
fn zero_iterations_free_path_f64() {
    zero_iterations_preserves_free_particle_path::<f64>();
}
