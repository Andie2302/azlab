//! Tests der Strömung: Quellenfreiheit, Ränder, Auftrieb und die
//! Teilchenflüssigkeit.

use super::{Edge, Edges, GridFluid, ParticleFluid};
use crate::geometry::Aabb;
use azmath::Vector;
use azmath::scalar::RealElementary;

fn v<T: RealElementary>(x: f64, y: f64) -> Vector<T, 2> {
    Vector([T::from_f64(x), T::from_f64(y)])
}

/// Toleranz, die zur Genauigkeit des Zahltyps passt.
fn loose<T: RealElementary>() -> T {
    T::EPSILON.sqrt()
}

fn rest<T: RealElementary>() {
    let f = T::from_f64;
    let mut fluid: GridFluid<T> = GridFluid::new(12, 9, f(0.1), Edges::closed());
    assert_eq!((fluid.nx(), fluid.ny()), (12, 9));
    fluid.add_solid_disc(v(0.6, 0.4), f(0.15));
    assert!(fluid.is_solid(5, 4) && !fluid.is_solid(0, 0));
    for _ in 0..5 {
        let stats = fluid.step(f(0.02));
        assert!(stats.converged);
    }
    // Ohne Antrieb bleibt alles in Ruhe.
    for j in 0..9 {
        for i in 0..12 {
            assert_eq!(fluid.velocity(i, j), v(0.0, 0.0));
            assert_eq!(fluid.smoke_at(i, j), T::ZERO);
        }
    }
    // Ungültige Schrittweite und Zellen außerhalb ändern nichts.
    fluid.add_smoke(99, 0, T::ONE);
    fluid.add_smoke(5, 4, T::ONE);
    fluid.step(T::ZERO);
    fluid.step(T::NAN);
    assert_eq!(fluid.smoke().iter().fold(T::ZERO, |sum, &value| sum + value), T::ZERO);
    assert_eq!(fluid.velocity(99, 0), v(0.0, 0.0));
}

fn divergence_free<T: RealElementary>() {
    let f = T::from_f64;
    let mut fluid: GridFluid<T> = GridFluid::new(24, 16, f(0.05), Edges::closed());
    fluid.params.pressure_tolerance = loose::<T>() * f(0.01);
    fluid.params.viscosity = f(0.01);
    fluid.add_solid_disc(v(0.6, 0.4), f(0.12));
    // Kräftig umrühren.
    for j in 0..16 {
        for i in 0..24 {
            let (x, y) = (i as f64 * 0.7, j as f64 * 1.3);
            fluid.add_velocity(i, j, v(x.sin() + 0.5 * y.cos(), (x + y).cos()));
        }
    }
    assert!(fluid.max_divergence() > f(1.0));
    for _ in 0..4 {
        let stats = fluid.step(f(0.01));
        assert!(stats.converged);
        assert!(fluid.max_divergence() < loose::<T>() * f(20.0));
    }
    // Im Hindernis und an seinen Kanten steht die Flüssigkeit.
    assert_eq!(fluid.velocity(12, 8), v(0.0, 0.0));
    assert_eq!(fluid.divergence(12, 8), T::ZERO);
    assert_eq!(fluid.vorticity(12, 8), T::ZERO);
    // Die Bewegung ist noch da und endlich.
    let energy = (0..24 * 16).fold(T::ZERO, |sum, index| sum + fluid.velocity(index % 24, index / 24).norm_sqr());
    assert!(energy > f(1.0) && energy.is_finite());
}

fn plume<T: RealElementary>() {
    let f = T::from_f64;
    let mut fluid: GridFluid<T> = GridFluid::new(20, 30, f(0.05), Edges::closed());
    fluid.params.buoyancy = f(4.0);
    for j in 2..6 {
        for i in 8..12 {
            fluid.add_smoke(i, j, T::ONE);
        }
    }
    let height = |fluid: &GridFluid<T>| {
        let total = fluid.smoke().iter().fold(T::ZERO, |sum, &value| sum + value);
        let moment = (0..20 * 30).fold(T::ZERO, |sum, index| sum + fluid.smoke()[index] * f((index / 20) as f64));
        moment / total
    };
    let start = height(&fluid);
    for _ in 0..60 {
        fluid.step(f(0.02));
    }
    // Der Rauch ist gestiegen, nie negativ und nicht mehr geworden als das
    // Maximum am Anfang.
    assert!(height(&fluid) > start + f(3.0));
    assert!(fluid.smoke().iter().all(|&value| value >= T::ZERO && value <= T::ONE));
    assert!(fluid.max_divergence() < f(0.05));
    // Mit Schwund verschwindet er.
    fluid.params.dissipation = f(50.0);
    for _ in 0..20 {
        fluid.step(f(0.02));
    }
    assert!(fluid.smoke().iter().all(|&value| value < f(1.0e-6)));
}

fn channel<T: RealElementary>() {
    let f = T::from_f64;
    let mut fluid: GridFluid<T> = GridFluid::new(30, 10, f(0.1), Edges::channel(f(2.0)));
    fluid.params.pressure_tolerance = loose::<T>() * f(0.01);
    for _ in 0..40 {
        fluid.step(f(0.02));
    }
    // Was links hereinströmt, strömt überall gleichmäßig durch.
    for j in 0..10 {
        for i in 0..30 {
            assert!((fluid.velocity(i, j) - v(2.0, 0.0)).norm() < loose::<T>() * f(50.0));
        }
    }
    // Mit einem Hindernis muss die Strömung ausweichen: Sie wird neben ihm
    // schneller und bleibt quellenfrei.
    fluid.add_solid_disc(v(1.0, 0.5), f(0.2));
    for _ in 0..40 {
        fluid.step(f(0.02));
    }
    assert!(fluid.velocity(10, 1)[0] > f(2.2));
    assert!(fluid.max_divergence() < loose::<T>() * f(50.0));
    // Offener Rand oben, Einlass unten: Es strömt senkrecht durch.
    let edges = Edges { left: Edge::Wall, right: Edge::Wall, bottom: Edge::Inflow(f(1.0)), top: Edge::Open };
    let mut column: GridFluid<T> = GridFluid::new(6, 12, f(0.1), edges);
    for _ in 0..20 {
        column.step(f(0.02));
    }
    assert!((column.velocity(3, 6) - v(0.0, 1.0)).norm() < f(1.0e-3));
}

fn tank<T: RealElementary>() -> ParticleFluid<T> {
    ParticleFluid::new(Aabb { min: v(0.0, 0.0), max: v(2.0, 1.2) }, T::from_f64(0.05))
}

fn centre<T: RealElementary>(fluid: &ParticleFluid<T>) -> Vector<T, 2> {
    let sum = fluid.positions().iter().fold(Vector::ZERO, |sum, &point| sum + point);
    sum * T::from_f64(fluid.positions().len() as f64).recip()
}

fn pool<T: RealElementary>() {
    let f = T::from_f64;
    let mut fluid: ParticleFluid<T> = tank();
    assert!(!fluid.walls().is_empty());
    // Becken bis zur Höhe 0,4 gefüllt: Das Wasser bleibt, wo es ist.
    assert_eq!(fluid.add_block(v(0.0, 0.0), v(2.0, 0.4)), 40 * 8);
    let start = centre(&fluid);
    for _ in 0..240 {
        fluid.step(f(1.0 / 60.0));
    }
    let moved = centre(&fluid) - start;
    assert!(moved.norm() < f(0.03), "Schwerpunkt um {:?} verschoben", moved.norm().to_f64());
    assert!(fluid.mean_compression() < f(0.05));
    // Nach dem Einschwingen zwei weitere Sekunden beobachten, nicht nur ein Bild.
    for _ in 0..120 {
        fluid.step(f(1.0 / 60.0));
        for (low, high) in [(0.0, 0.15), (0.6, 1.4), (1.85, 2.0)] {
            let mut count = 0usize;
            let (mut sum, mut maximum) = (T::ZERO, T::ZERO);
            for (point, velocity) in fluid.positions().iter().zip(fluid.velocities()) {
                if point[0] >= f(low) && point[0] < f(high) {
                    let speed = velocity.norm();
                    assert!(speed.is_finite());
                    count += 1;
                    sum += speed;
                    maximum = maximum.max(speed);
                }
            }
            assert!(count > 0, "leere Messregion {low}..{high}");
            assert!(sum / f(count as f64) < f(0.02), "Randströmung {low}..{high}");
            assert!(maximum < f(0.1), "schnelles Einzelteilchen {low}..{high}: {:?}", maximum.to_f64());
        }
    }
    // Ungültige Schrittweite ändert nichts.
    let before = fluid.positions()[0];
    fluid.step(T::ZERO);
    fluid.step(T::NAN);
    assert_eq!(fluid.positions()[0], before);
}

fn dam_break<T: RealElementary>() {
    let f = T::from_f64;
    let mut fluid: ParticleFluid<T> = tank();
    fluid.add_obstacle(v(1.2, 0.0), v(1.3, 0.2));
    assert_eq!(fluid.obstacles().len(), 1);
    // Wassersäule links; nichts landet im Hindernis oder außerhalb.
    let count = fluid.add_block(v(0.0, 0.0), v(0.5, 1.0));
    assert_eq!(count, 10 * 20);
    assert_eq!(fluid.add_block(v(1.2, 0.0), v(1.3, 0.2)), 0);
    assert_eq!(fluid.add_block(v(5.0, 5.0), v(6.0, 6.0)), 0);
    let start = centre(&fluid);
    let bounds = fluid.bounds();
    for _ in 0..150 {
        fluid.step(f(1.0 / 60.0));
        // Kein Teilchen verlässt das Becken oder dringt ins Hindernis ein.
        for &point in fluid.positions() {
            assert!(bounds.contains(point));
            assert!(!fluid.obstacles()[0].contains(point));
        }
    }
    // Die Säule ist zusammengesackt und nach rechts geflossen.
    let end = centre(&fluid);
    assert!(end[1] < start[1] - f(0.2));
    assert!(end[0] > start[0] + f(0.2));
    // Die Flüssigkeit lässt sich kaum zusammendrücken.
    assert!(fluid.mean_compression() < f(0.08), "Verdichtung {:?}", fluid.mean_compression().to_f64());
    assert!(fluid.densities().iter().all(|density| density.is_finite()));
    assert_eq!(fluid.spacing(), f(0.05));
}

#[test]
fn grid_f32() {
    rest::<f32>();
    divergence_free::<f32>();
    plume::<f32>();
    channel::<f32>();
}

#[test]
fn grid_f64() {
    rest::<f64>();
    divergence_free::<f64>();
    plume::<f64>();
    channel::<f64>();
}

#[test]
fn particles_f32() {
    pool::<f32>();
    dam_break::<f32>();
}

#[test]
fn particles_f64() {
    pool::<f64>();
    dam_break::<f64>();
}

#[cfg(feature = "std")]
#[test]
#[ignore = "Parametervergleich über 20 simulierte Sekunden; mit --ignored --nocapture starten"]
fn particle_resolution_and_obstacle_sweep() {
    // Messung statt behaupteter universeller Ruheschwelle: dünne Wände sind
    // ausdrücklich eingeschränkt. Sicherheit gilt bei jeder getesteten Dicke.
    for spacing in [0.05f64, 0.1] {
        for layers in [1.0, 3.0, 5.0] {
            let mut fluid = ParticleFluid::new(
                Aabb { min: v(0.0, 0.0), max: v(2.0, 1.2) }, spacing);
            fluid.add_obstacle(v(1.0, 0.0), v(1.0 + layers * spacing, 0.2));
            fluid.add_block(v(0.0, 0.0), v(2.0, 0.4));
            let mut peak = 0.0f64;
            for frame in 0..1200 {
                fluid.step(1.0 / 60.0);
                for (&p, velocity) in fluid.positions().iter().zip(fluid.velocities()) {
                    assert!(fluid.bounds().contains(p));
                    assert!(!fluid.obstacles()[0].contains(p));
                    let speed = velocity.norm();
                    assert!(speed.is_finite());
                    if frame >= 900 { peak = peak.max(speed); }
                }
            }
            std::println!("spacing={spacing}, obstacle_layers={layers}, maximum_15_to_20s={peak} m/s");
        }
    }
}

#[test]
#[ignore = "Kleines volles Kontrollbecken, 20 simulierte Sekunden"]
fn small_full_pool_defaults() {
    let mut fluid = ParticleFluid::<f32>::new(
        Aabb { min: v(0.0, 0.0), max: v(1.0, 1.0) }, 0.1);
    assert_eq!(fluid.params.substeps, 2);
    assert!((fluid.params.anti_clumping - 0.002).abs() < 1.0e-7);
    assert_eq!(fluid.add_block(v(0.0, 0.0), v(1.0, 1.0)), 100);
    for frame in 0..1200 {
        fluid.step(1.0 / 60.0);
        if frame < 900 { continue; }
        for (low, high) in [(0.0, 0.15), (0.3, 0.7), (0.85, 1.0)] {
            let speeds: alloc::vec::Vec<_> = fluid.positions().iter().zip(fluid.velocities())
                .filter(|(p, _)| p[0] >= low && p[0] < high)
                .map(|(_, v)| v.norm()).collect();
            assert!(!speeds.is_empty());
            assert!(speeds.iter().all(|s| s.is_finite()));
            let mean = speeds.iter().sum::<f32>() / speeds.len() as f32;
            let max = speeds.iter().copied().fold(0.0f32, f32::max);
            assert!(mean < 0.02, "frame={frame}, region={low}..{high}, mean={mean}");
            assert!(max < 0.1, "frame={frame}, region={low}..{high}, max={max}");
        }
    }
}

#[test]
fn fluid_motion_restore_is_atomic_and_can_continue() {
    let mut fluid = tank::<f64>();
    fluid.add_block(v(0.0, 0.0), v(2.0, 0.4));
    let positions = fluid.positions().to_vec();
    let velocities = fluid.velocities().to_vec();
    let mut bad = positions.clone();
    bad[0] = v(-1.0, 0.0);
    assert_eq!(fluid.restore_motion(&bad, &velocities), Err(crate::MotionError::InvalidPosition));
    assert_eq!(fluid.positions(), positions);
    let mut bad_velocity = velocities.clone();
    bad_velocity[0][0] = f64::INFINITY;
    assert_eq!(fluid.restore_motion(&positions, &bad_velocity), Err(crate::MotionError::NonFinite));
    fluid.step(1.0 / 60.0);
    assert!(fluid.restore_motion(&positions, &velocities).is_ok());
    let mut expected = tank::<f64>();
    expected.add_block(v(0.0, 0.0), v(2.0, 0.4));
    for _ in 0..3 {
        fluid.step(1.0 / 60.0);
        expected.step(1.0 / 60.0);
    }
    assert_eq!(fluid.positions(), expected.positions());
    assert_eq!(fluid.velocities(), expected.velocities());
}
