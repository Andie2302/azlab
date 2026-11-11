//! Tests der Gitteroperatoren an Standardproblemen mit bekannter Lösung.

#![cfg(any(feature = "std", feature = "libm"))]

use super::{
    Boundary, Derivative, Laplacian, LinearSolver, Scheme, UniformGrid, poisson_multigrid_1d,
    poisson_multigrid_2d, theta_step,
};
use crate::linalg::{AffineOperator, LinalgError, LinearOperator};
use crate::ode::{euler_slice, rk4_slice};
use crate::scalar::RealElementary;
use crate::solve::{SolveError, Tolerance, conjugate_gradient};
use crate::testutil::close;

/// Größter Betrag der Differenz zweier Felder.
fn max_difference<T: RealElementary>(a: &[T], b: &[T]) -> T {
    a.iter().zip(b).fold(T::ZERO, |worst, (&x, &y)| worst.max((x - y).abs()))
}

/// Inneres Gitter auf (0, 1) mit `N` Punkten und Rand bei 0 und 1.
fn interior<T: RealElementary, const N: usize>() -> Option<(UniformGrid<T, 1>, T)> {
    let h = T::ONE / T::from_f64((N + 1) as f64);
    UniformGrid::new([N], [h], [h]).map(|grid| (grid, h))
}

fn grids<T: RealElementary>() {
    let f = T::from_f64;
    let grid = UniformGrid::new([3, 4], [f(1.0), f(-1.0)], [f(0.5), f(0.25)]);
    assert!(grid.is_some_and(|g| {
        g.len() == 12
            && !g.is_empty()
            && g.shape() == [3, 4]
            && g.spacing() == [f(0.5), f(0.25)]
            && g.flat([2, 3]) == Some(11)
            && g.flat([1, 2]) == Some(6)
            && g.flat([3, 0]).is_none()
            && g.flat([0, 4]).is_none()
            && g.index_of(6) == Some([1, 2])
            && g.index_of(12).is_none()
            && g.point([2, 1]) == [f(2.0), f(-0.75)]
            && (0..12).all(|flat| g.index_of(flat).and_then(|index| g.flat(index)) == Some(flat))
    }));
    assert!(UniformGrid::new([3], [T::ZERO], [T::ZERO]).is_none());
    assert!(UniformGrid::new([3], [T::ZERO], [T::NAN]).is_none());
    assert!(UniformGrid::<T, 2>::new([usize::MAX, 2], [T::ZERO; 2], [T::ONE; 2]).is_none());
    assert!(UniformGrid::<T, 1>::new([0], [T::ZERO], [T::ONE]).is_some_and(|g| g.is_empty()));

    // Abtasten einer Funktion.
    let mut field = [T::ZERO; 12];
    assert!(grid.is_some_and(|g| g.sample(&mut field, |[x, y]| x + y)));
    assert_eq!(field[6], f(1.5) + f(-0.5));
    let mut wrong = [T::ZERO; 5];
    assert!(grid.is_some_and(|g| !g.sample(&mut wrong, |_| T::ONE)));
}

fn stencils<T: RealElementary>() {
    let f = T::from_f64;
    let pi = T::PI;

    // Δ sin(πx) = −π² sin(πx); der Fehler fällt mit h² (Faktor 4 je Halbierung).
    let error = |points: usize| -> Option<T> {
        let h = T::ONE / f((points + 1) as f64);
        let grid = UniformGrid::new([points], [h], [h])?;
        let mut u = [T::ZERO; 63];
        let mut image = [T::ZERO; 63];
        let (u, image) = (u.get_mut(..points)?, image.get_mut(..points)?);
        grid.sample(u, |[x]| (pi * x).sin());
        Laplacian::new(&grid, [Boundary::Dirichlet]).apply(u, image).ok()?;
        Some(u.iter().zip(image.iter()).fold(T::ZERO, |worst, (&value, &second)| {
            worst.max((second + pi * pi * value).abs())
        }))
    };
    assert!(match (error(15), error(31)) {
        (Some(coarse), Some(fine)) => fine < coarse * f(0.3) && fine < f(0.05),
        _ => false,
    });

    // 2D: Δ sin(πx)·sin(πy) = −2π² sin(πx)·sin(πy).
    let h = T::ONE / f(16.0);
    let grid = UniformGrid::new([15, 15], [h, h], [h, h]);
    let mut u = [T::ZERO; 225];
    let mut image = [T::ZERO; 225];
    assert!(grid.is_some_and(|g| {
        let laplace = Laplacian::new(&g, [Boundary::Dirichlet; 2]);
        g.sample(&mut u, |[x, y]| (pi * x).sin() * (pi * y).sin())
            && laplace.apply(&u, &mut image).is_ok()
            && laplace.rows() == 225
            && laplace.cols() == 225
    }));
    let worst = u.iter().zip(&image).fold(T::ZERO, |w, (&value, &second)| {
        w.max((second + f(2.0) * pi * pi * value).abs())
    });
    assert!(worst < f(0.1));

    // Symmetrie mit Dirichlet-Rand: ⟨Δu, v⟩ = ⟨u, Δv⟩.
    let mut v = [T::ZERO; 225];
    let mut image_v = [T::ZERO; 225];
    assert!(grid.is_some_and(|g| {
        g.sample(&mut v, |[x, y]| x * x - y + f(0.3))
            && Laplacian::new(&g, [Boundary::Dirichlet; 2]).apply(&v, &mut image_v).is_ok()
    }));
    let dot = |a: &[T], b: &[T]| a.iter().zip(b).fold(T::ZERO, |sum, (&x, &y)| sum + x * y);
    assert!(close(dot(&image, &v), dot(&u, &image_v), 4096.0));

    // Neumann und periodisch: Konstanten liegen im Kern, die Summe von Δu ist null.
    let line = UniformGrid::new([8], [T::ZERO], [f(0.125)]);
    for boundary in [Boundary::Neumann, Boundary::Periodic] {
        let ones = [T::ONE; 8];
        let mut out = [T::ONE; 8];
        let mut wave = [T::ZERO; 8];
        let mut second = [T::ZERO; 8];
        assert!(line.is_some_and(|g| {
            let laplace = Laplacian::new(&g, [boundary]);
            g.sample(&mut wave, |[x]| (f(2.0) * pi * x).sin() + x * x)
                && laplace.apply(&ones, &mut out).is_ok()
                && laplace.apply(&wave, &mut second).is_ok()
        }));
        assert!(out.iter().all(|&value| value == T::ZERO));
        assert!(second.iter().fold(T::ZERO, |sum, &value| sum + value).abs() < f(1.0e-3));
    }

    // Erste Ableitung, periodisch: (sin 2πx)' = 2π cos 2πx; die zentrale
    // Formel ist genauer als die einseitigen, alle erhalten die Masse.
    let periodic = UniformGrid::new([64], [T::ZERO], [T::ONE / f(64.0)]);
    let mut wave = [T::ZERO; 64];
    let mut exact = [T::ZERO; 64];
    assert!(periodic.is_some_and(|g| {
        g.sample(&mut wave, |[x]| (f(2.0) * pi * x).sin())
            && g.sample(&mut exact, |[x]| f(2.0) * pi * (f(2.0) * pi * x).cos())
    }));
    let deviation = |scheme: Scheme| -> Option<(T, T)> {
        let mut slope = [T::ZERO; 64];
        Derivative::new(&periodic?, 0, scheme, Boundary::Periodic)?.apply(&wave, &mut slope).ok()?;
        Some((max_difference(&slope, &exact), slope.iter().fold(T::ZERO, |sum, &s| sum + s)))
    };
    assert!(match (deviation(Scheme::Central), deviation(Scheme::Forward), deviation(Scheme::Backward)) {
        (Some((central, mass_c)), Some((forward, mass_f)), Some((backward, mass_b))) => {
            central < f(0.02)
                && central < forward * f(0.2)
                && central < backward * f(0.2)
                && mass_c.abs() < f(1.0e-3)
                && mass_f.abs() < f(1.0e-3)
                && mass_b.abs() < f(1.0e-3)
        }
        _ => false,
    });
    // Ableitung in einer Richtung eines 2D-Felds lässt die andere unberührt.
    let plane = UniformGrid::new([4, 5], [T::ZERO; 2], [f(0.5), f(0.25)]);
    let mut field = [T::ZERO; 20];
    let mut slope = [T::ZERO; 20];
    assert!(plane.is_some_and(|g| {
        g.sample(&mut field, |[x, y]| f(3.0) * x + f(7.0) * y)
            && Derivative::new(&g, 0, Scheme::Central, Boundary::Neumann)
                .is_some_and(|d| d.apply(&field, &mut slope).is_ok() && d.rows() == 20 && d.cols() == 20)
            && Derivative::new(&g, 2, Scheme::Central, Boundary::Neumann).is_none()
    }));
    // Im Inneren (Zeilen 1 und 2) ist ∂/∂x exakt 3.
    assert!(slope[5..15].iter().all(|&value| close(value, f(3.0), 16.0)));

    // Falsche Längen.
    let mut short = [T::ZERO; 3];
    assert!(line.is_some_and(|g| {
        Laplacian::new(&g, [Boundary::Dirichlet]).apply(&[T::ONE; 8], &mut short)
            == Err(LinalgError::DimensionMismatch)
    }));
}

fn poisson<T: RealElementary>() {
    let f = T::from_f64;
    let pi = T::PI;
    let tol = Tolerance::<T>::default().with_rel_x(T::EPSILON.sqrt()).with_max_iterations(500);

    // −u'' = π² sin(πx), u(0) = u(1) = 0 → u = sin(πx).
    let Some((grid, h)) = interior::<T, 31>() else {
        return assert!(interior::<T, 31>().is_some());
    };
    let mut rhs = [T::ZERO; 31];
    let mut exact = [T::ZERO; 31];
    grid.sample(&mut rhs, |[x]| pi * pi * (pi * x).sin());
    grid.sample(&mut exact, |[x]| (pi * x).sin());

    // Konjugierte Gradienten auf −Δ.
    let laplace = Laplacian::new(&grid, [Boundary::Dirichlet]);
    let negative = AffineOperator::new(&laplace, T::ZERO, -T::ONE);
    let mut u = [T::ZERO; 31];
    let mut work = [T::ZERO; 155];
    assert!(conjugate_gradient(&negative, &rhs, &mut u, &mut work, &tol).is_ok());
    assert!(max_difference(&u, &exact) < f(2.0e-3));

    // Mehrgitter: dieselbe Lösung in wenigen Zyklen.
    let mut v = [T::ZERO; 31];
    assert!(poisson_multigrid_1d(&rhs, &mut v, h, &mut work, &tol).is_ok_and(|s| s.iterations <= 12));
    assert!(max_difference(&v, &u) < f(1.0e-3));
    // Fehlerfälle.
    assert_eq!(
        poisson_multigrid_1d(&rhs[..30], &mut v[..30], h, &mut work, &tol),
        Err(SolveError::InvalidInput)
    );
    assert_eq!(
        poisson_multigrid_1d(&rhs, &mut v, h, &mut work[..100], &tol),
        Err(SolveError::InvalidInput)
    );
    assert!(matches!(
        poisson_multigrid_1d(&rhs, &mut [T::ZERO; 31], h, &mut work, &tol.with_max_iterations(0)),
        Err(SolveError::NotConverged(_))
    ));

    // 2D: −Δu = 2π² sin(πx)·sin(πy) → u = sin(πx)·sin(πy).
    let h = T::ONE / f(16.0);
    let plane = UniformGrid::new([15, 15], [h, h], [h, h]);
    let mut rhs = [T::ZERO; 225];
    let mut exact = [T::ZERO; 225];
    assert!(plane.is_some_and(|g| {
        g.sample(&mut rhs, |[x, y]| f(2.0) * pi * pi * (pi * x).sin() * (pi * y).sin())
            && g.sample(&mut exact, |[x, y]| (pi * x).sin() * (pi * y).sin())
    }));
    let mut u = [T::ZERO; 225];
    let mut work = [T::ZERO; 1125];
    assert!(poisson_multigrid_2d(&rhs, &mut u, 15, h, &mut work, &tol).is_ok_and(|s| s.iterations <= 15));
    assert!(max_difference(&u, &exact) < f(5.0e-3));
    // Gegenprobe mit konjugierten Gradienten.
    let mut v = [T::ZERO; 225];
    assert!(plane.is_some_and(|g| {
        let laplace = Laplacian::new(&g, [Boundary::Dirichlet; 2]);
        let negative = AffineOperator::new(&laplace, T::ZERO, -T::ONE);
        conjugate_gradient(&negative, &rhs, &mut v, &mut work, &tol).is_ok()
    }));
    assert!(max_difference(&u, &v) < f(1.0e-3));
    assert_eq!(
        poisson_multigrid_2d(&rhs, &mut u, 14, h, &mut work, &tol),
        Err(SolveError::InvalidInput)
    );
}

fn evolution<T: RealElementary>() {
    let f = T::from_f64;
    let pi = T::PI;
    let tol = Tolerance::<T>::default().with_rel_x(T::EPSILON.sqrt()).with_max_iterations(500);
    let Some((grid, _)) = interior::<T, 31>() else {
        return assert!(interior::<T, 31>().is_some());
    };
    let laplace = Laplacian::new(&grid, [Boundary::Dirichlet]);
    let mut start = [T::ZERO; 31];
    grid.sample(&mut start, |[x]| (pi * x).sin());

    // Wärmeleitung u_t = u_xx: Die Grundschwingung klingt wie e^(−π²t) ab.
    let (steps, dt) = (50, f(0.002));
    let decay = (-pi * pi * dt * f(f64::from(steps))).exp();
    for (theta, limit) in [(f(0.5), f(2.0e-3)), (T::ONE, f(2.0e-2))] {
        let mut u = start;
        let mut work = [T::ZERO; 217];
        for _ in 0..steps {
            assert!(theta_step(&laplace, &mut u, dt, theta, LinearSolver::ConjugateGradient, &mut work, &tol).is_ok());
        }
        let expected = start.map(|value| value * decay);
        assert!(max_difference(&u, &expected) < limit);
    }
    // BiCGSTAB liefert denselben Schritt wie konjugierte Gradienten.
    let (mut a, mut b) = (start, start);
    let mut work = [T::ZERO; 217];
    assert!(theta_step(&laplace, &mut a, dt, f(0.5), LinearSolver::ConjugateGradient, &mut work, &tol).is_ok());
    assert!(theta_step(&laplace, &mut b, dt, f(0.5), LinearSolver::BiCgStab, &mut work, &tol).is_ok());
    assert!(max_difference(&a, &b) < f(1.0e-3));
    assert_eq!(
        theta_step(&laplace, &mut a, dt, f(0.5), LinearSolver::BiCgStab, &mut work[..100], &tol),
        Err(SolveError::InvalidInput)
    );

    // Wellengleichung u_tt = u_xx als System (u, v): Die Grundschwingung
    // kehrt nach der Periode 2 zum Anfangszustand zurück.
    let mut state = [T::ZERO; 62];
    state[..31].copy_from_slice(&start);
    let mut work = [T::ZERO; 186];
    let wave = |_t: T, y: &[T], dy: &mut [T]| {
        let (u, v) = y.split_at(31);
        let (du, dv) = dy.split_at_mut(31);
        du.copy_from_slice(v);
        // Mit der diskreten Schwingungsdauer der Grundschwingung vergleichen
        // wäre genauer; für den Test genügt die kontinuierliche.
        let _ = laplace.apply(u, dv);
    };
    let (steps, dt) = (400_u32, f(0.005));
    for step in 0..steps {
        assert!(rk4_slice(wave, dt * f(f64::from(step)), &mut state, dt, &mut work).is_ok());
    }
    assert!(max_difference(&state[..31], &start) < f(2.0e-2));

    // Transport u_t + u_x = 0, periodisch, Aufwind-Formel mit Euler-Schritten:
    // Die Masse bleibt erhalten, es entstehen keine neuen Extrema.
    let ring = UniformGrid::new([64], [T::ZERO], [T::ONE / f(64.0)]);
    let mut u = [T::ZERO; 64];
    assert!(ring.is_some_and(|g| g.sample(&mut u, |[x]| (f(2.0) * pi * x).sin() + f(2.0))));
    let mass = u.iter().fold(T::ZERO, |sum, &value| sum + value);
    let upwind = ring.and_then(|g| Derivative::new(&g, 0, Scheme::Backward, Boundary::Periodic));
    let mut work = [T::ZERO; 64];
    let dt = f(0.5) / f(64.0);
    assert!(upwind.is_some_and(|derivative| {
        (0..100_u32).all(|_| {
            euler_slice(
                |_t: T, y: &[T], dy: &mut [T]| {
                    let _ = derivative.apply(y, dy);
                    dy.iter_mut().for_each(|slope| *slope = -*slope);
                },
                T::ZERO,
                &mut u,
                dt,
                &mut work,
            )
            .is_ok()
        })
    }));
    assert!(close(u.iter().fold(T::ZERO, |sum, &value| sum + value), mass, 4096.0));
    assert!(u.iter().all(|&value| value >= T::ONE - f(1.0e-4) && value <= f(3.0) + f(1.0e-4)));
    // Zu kleiner Arbeitsspeicher.
    assert_eq!(
        rk4_slice(|_t: T, _y: &[T], _dy: &mut [T]| {}, T::ZERO, &mut u, dt, &mut work),
        Err(LinalgError::DimensionMismatch)
    );
    assert_eq!(
        euler_slice(|_t: T, _y: &[T], _dy: &mut [T]| {}, T::ZERO, &mut u, dt, &mut work[..10]),
        Err(LinalgError::DimensionMismatch)
    );
}

#[test]
fn pde_f32() {
    grids::<f32>();
    stencils::<f32>();
    poisson::<f32>();
    evolution::<f32>();
}

#[test]
fn pde_f64() {
    grids::<f64>();
    stencils::<f64>();
    poisson::<f64>();
    evolution::<f64>();
}
