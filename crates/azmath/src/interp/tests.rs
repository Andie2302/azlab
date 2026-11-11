//! Tests der Interpolation.

use super::{Barycentric, CubicSpline, InterpError, lerp, linear, natural_spline, spline_eval};
use crate::scalar::Real;
use crate::testutil::close;

fn basics<T: Real>() {
    let f = T::from_f64;
    let xs = [f(0.0), f(1.0), f(3.0), f(4.0)];
    let ys = [f(1.0), f(3.0), f(-1.0), f(0.0)];

    assert_eq!(lerp(f(2.0), f(6.0), f(0.25)), f(3.0));
    // Linear: trifft die Punkte, liegt dazwischen auf der Strecke.
    assert_eq!(linear(&xs, &ys, f(1.0)), Some(f(3.0)));
    assert_eq!(linear(&xs, &ys, f(0.0)), Some(f(1.0)));
    assert_eq!(linear(&xs, &ys, f(4.0)), Some(f(0.0)));
    assert_eq!(linear(&xs, &ys, f(2.0)), Some(f(1.0)));
    assert_eq!(linear(&xs, &ys, f(0.5)), Some(f(2.0)));
    assert_eq!(linear(&xs, &ys, f(-0.1)), None);
    assert_eq!(linear(&xs, &ys, f(4.1)), None);
    assert_eq!(linear(&xs, &ys, T::NAN), None);
    assert_eq!(linear(&xs, &ys[..3], f(1.0)), None);
    assert_eq!(linear(&xs[..1], &ys[..1], f(0.0)), None);

    // Spline: trifft die Punkte, ist an den Enden „natürlich“, und stimmt
    // für eine Gerade mit ihr überein.
    let spline = CubicSpline::natural(xs, ys);
    assert!(spline.is_ok_and(|s| {
        xs.iter().zip(ys).all(|(&x, y)| s.eval(x).is_some_and(|v| close(v, y, 16.0)))
            && s.eval(f(-0.1)).is_none()
            && s.derivative(f(5.0)).is_none()
            && s.eval(f(2.0)).is_some()
    }));
    let line = CubicSpline::natural(xs, xs.map(|x| f(2.0) * x + T::ONE));
    assert!(line.is_ok_and(|s| {
        s.eval(f(2.5)).is_some_and(|v| close(v, f(6.0), 16.0))
            && s.derivative(f(2.5)).is_some_and(|d| close(d, f(2.0), 64.0))
    }));
    // Ableitung passt zum Differenzenquotienten des Splines.
    assert!(spline.is_ok_and(|s| {
        let h = f(1.0e-2);
        match (s.eval(f(2.0) + h), s.eval(f(2.0) - h), s.derivative(f(2.0))) {
            (Some(right), Some(left), Some(slope)) => {
                ((right - left) / (h + h) - slope).abs() < f(1.0e-2)
            }
            _ => false,
        }
    }));
    // Slice-Fassung liefert dasselbe.
    let mut m = [T::ZERO; 4];
    let mut work = [T::ZERO; 4];
    assert!(natural_spline(&xs, &ys, &mut m, &mut work).is_ok());
    assert_eq!((m[0], m[3]), (T::ZERO, T::ZERO));
    assert_eq!(spline_eval(&xs, &ys, &m, f(2.0)), spline.ok().and_then(|s| s.eval(f(2.0))));
    // Fehlerfälle.
    assert_eq!(
        CubicSpline::natural([f(0.0), f(0.0)], [f(1.0), f(2.0)]).err(),
        Some(InterpError::NotIncreasing)
    );
    assert_eq!(CubicSpline::natural([f(0.0)], [f(1.0)]).err(), Some(InterpError::TooFewPoints));
    assert_eq!(
        natural_spline(&xs, &ys, &mut m[..3], &mut work),
        Err(InterpError::LengthMismatch)
    );
    // Zwei Punkte: der Spline ist die Verbindungsstrecke.
    assert!(CubicSpline::natural([f(0.0), f(2.0)], [f(1.0), f(5.0)])
        .is_ok_and(|s| s.eval(f(0.5)).is_some_and(|v| close(v, f(2.0), 8.0))));

    // Polynom durch drei Punkte der Parabel x² ist die Parabel.
    let parabola = Barycentric::new([f(-1.0), f(0.5), f(2.0)], [f(1.0), f(0.25), f(4.0)]);
    assert!(parabola.is_ok_and(|p| {
        close(p.eval(f(3.0)), f(9.0), 64.0)
            && close(p.eval(f(1.0)), f(1.0), 16.0)
            && p.eval(f(0.5)) == f(0.25)
    }));
    assert_eq!(
        Barycentric::new([f(1.0), f(1.0)], [f(0.0), f(1.0)]).err(),
        Some(InterpError::NotIncreasing)
    );
    assert_eq!(Barycentric::<T, 0>::new([], []).err(), Some(InterpError::TooFewPoints));
    assert!(Barycentric::new([f(1.0)], [f(7.0)]).is_ok_and(|p| p.eval(f(3.0)) == f(7.0)));
}

#[test]
fn basics_f32() {
    basics::<f32>();
}

#[test]
fn basics_f64() {
    basics::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod chebyshev {
    use super::*;
    use crate::interp::{ChebyshevSeries, chebyshev_nodes};
    use crate::scalar::RealElementary;

    fn check<T: RealElementary>() {
        let f = T::from_f64;
        // Stützstellen: aufsteigend, im Intervall, symmetrisch zur Mitte.
        let nodes: [T; 5] = chebyshev_nodes(f(1.0), f(3.0));
        assert!(nodes.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(nodes[0] > f(1.0) && nodes[4] < f(3.0));
        assert!(close(nodes[2], f(2.0), 8.0) && close(nodes[0] + nodes[4], f(4.0), 8.0));

        // Ein Polynom vom Grad < N wird exakt wiedergegeben.
        let cubic = |x: T| x * x * x - f(2.0) * x + T::ONE;
        let series = ChebyshevSeries::<T, 6>::fit(cubic, f(-1.0), f(2.0));
        for x in [f(-1.0), f(-0.3), f(0.7), f(2.0)] {
            assert!(close(series.eval(x), cubic(x), 64.0));
        }
        // Auf [−1, 1] ist x² = (T₀ + T₂)/2.
        let square = ChebyshevSeries::<T, 4>::fit(|x| x * x, -T::ONE, T::ONE);
        let c = square.coefficients();
        assert!(close(c[0], f(0.5), 16.0) && close(c[2], f(0.5), 16.0));
        assert!(close(c[1], T::ZERO, 16.0) && close(c[3], T::ZERO, 16.0));

        // Glatte Funktion: wenige Glieder reichen für hohe Genauigkeit.
        let smooth = ChebyshevSeries::<T, 16>::fit(|x| x.exp(), T::ZERO, f(2.0));
        let limit = T::EPSILON.sqrt();
        for x in [f(0.0), f(0.4), f(1.3), f(2.0)] {
            assert!((smooth.eval(x) - x.exp()).abs() <= limit);
        }
        // Polynom an den Tschebyschow-Stützstellen interpolieren.
        let nodes: [T; 8] = chebyshev_nodes(T::ZERO, f(2.0));
        let poly = Barycentric::new(nodes, nodes.map(|x| x.sin()));
        assert!(poly.is_ok_and(|p| (p.eval(f(1.1)) - f(1.1).sin()).abs() <= f(1.0e-5)));
    }

    #[test]
    fn chebyshev_f32() {
        check::<f32>();
    }

    #[test]
    fn chebyshev_f64() {
        check::<f64>();
    }
}
