//! Tests der Zerlegungen. Erwartungswerte sind von Hand herleitbar; Löser
//! werden zusätzlich über Rückeinsetzen geprüft.

use crate::linalg::{LinalgError, Matrix, Vector};
use crate::number::{Complex, Dual};
use crate::scalar::{Algebra, Real};
use crate::testutil::close;

fn m<T: Real, const R: usize, const C: usize>(rows: [[f64; C]; R]) -> Matrix<T, R, C> {
    Matrix(rows.map(|row| row.map(T::from_f64)))
}

fn v<T: Real, const N: usize>(data: [f64; N]) -> Vector<T, N> {
    Vector(data.map(T::from_f64))
}

fn close_v<T: Real, const N: usize>(a: Vector<T, N>, b: Vector<T, N>, ulps: f64) -> bool {
    a.0.into_iter().zip(b.0).all(|(x, y)| close(x, y, ulps))
}

fn close_m<T: Real, const R: usize, const C: usize>(
    a: Matrix<T, R, C>,
    b: Matrix<T, R, C>,
    ulps: f64,
) -> bool {
    a.0.into_iter().zip(b.0).all(|(x, y)| close_v(Vector(x), Vector(y), ulps))
}

fn lu<T: Real>() {
    let id2 = Matrix::<T, 2, 2>::identity();

    // 2x + y = 3, x + 3y = 5  →  x = 0.8, y = 1.4
    let a = m::<T, 2, 2>([[2.0, 1.0], [1.0, 3.0]]);
    let b = v::<T, 2>([3.0, 5.0]);
    assert!(a.solve(b).is_ok_and(|x| close_v(x, v([0.8, 1.4]), 8.0) && close_v(a * x, b, 8.0)));
    assert_eq!(a.det(), T::from_f64(5.0));
    assert!(a.inverse().is_ok_and(|inv| close_m(inv * a, id2, 8.0) && close_m(a * inv, id2, 8.0)));
    assert!(a.lu().is_ok_and(|lu| close(lu.det(), T::from_f64(5.0), 8.0)
        && close_m(lu.solve_matrix(a), id2, 8.0)
        && close_m(lu.l() * lu.u(), a, 8.0)
        && lu.permutation() == [0, 1]));

    // Zeilentausch nötig.
    let swap = m::<T, 2, 2>([[0.0, 1.0], [1.0, 0.0]]);
    assert_eq!(swap.solve(v([2.0, 3.0])), Ok(v([3.0, 2.0])));
    assert!(swap.lu().is_ok_and(|lu| lu.permutation() == [1, 0] && lu.det() == -T::ONE));

    // Determinanten: 3×3 geschlossen und über LU, 4×4 mit Vorzeichenwechsel.
    let c = m::<T, 3, 3>([[2.0, 0.0, 1.0], [1.0, 3.0, 2.0], [1.0, 1.0, 4.0]]);
    assert_eq!(c.det(), T::from_f64(18.0));
    assert!(c.lu().is_ok_and(|lu| close(lu.det(), T::from_f64(18.0), 8.0)));
    let d = m::<T, 4, 4>([
        [0.0, 2.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 3.0, 0.0],
        [0.0, 0.0, 0.0, 4.0],
    ]);
    assert!(close(d.det(), T::from_f64(-24.0), 8.0));
    assert_eq!(Matrix::<T, 0, 0>::ZERO.det(), T::ONE);
    assert_eq!(m::<T, 1, 1>([[7.0]]).det(), T::from_f64(7.0));

    // P·A = L·U und Residuum für ein 3×3-System.
    let rhs = v::<T, 3>([1.0, 2.0, 3.0]);
    assert!(c.lu().is_ok_and(|lu| {
        let p = lu.permutation();
        let permuted = Matrix::<T, 3, 3>::from_fn(|i, j| c[(p[i], j)]);
        close_m(lu.l() * lu.u(), permuted, 8.0) && close_v(c * lu.solve(rhs), rhs, 16.0)
    }));

    // Singulär.
    let s = m::<T, 2, 2>([[1.0, 2.0], [2.0, 4.0]]);
    assert!(matches!(s.lu(), Err(LinalgError::Singular)));
    assert_eq!(s.det(), T::ZERO);
    assert_eq!(s.inverse(), Err(LinalgError::Singular));
    assert_eq!(s.solve(v([1.0, 1.0])), Err(LinalgError::Singular));
    assert_eq!(Matrix::<T, 4, 4>::ZERO.det(), T::ZERO);

    // Komplex: det [[i, 1], [1, i]] = i² − 1 = −2, x = A⁻¹·(1, 0) = (−i/2, 1/2).
    let i = Complex::<T>::I;
    let one = Complex::<T>::ONE;
    let z = Matrix([[i, one], [one, i]]);
    assert_eq!(z.det(), Complex::from_real(T::from_f64(-2.0)));
    assert!(z.solve(Vector([one, Complex::ZERO])).is_ok_and(|x| {
        close(x[0].re, T::ZERO, 8.0)
            && close(x[0].im, T::from_f64(-0.5), 8.0)
            && close(x[1].re, T::from_f64(0.5), 8.0)
            && close(x[1].im, T::ZERO, 8.0)
    }));

    // Ableitung durch das Lösen: A(p) = [[p, 1], [1, 2]], b = (1, 0).
    // x₀ = 2/(2p − 1), x₁ = −1/(2p − 1); bei p = 1: x = (2, −1),
    // dx/dp = (−4, 2).
    let k = |value: f64| Dual::constant(T::from_f64(value));
    let p = Dual::variable(T::ONE);
    let dual = Matrix([[p, k(1.0)], [k(1.0), k(2.0)]]);
    assert!(dual.solve(Vector([k(1.0), k(0.0)])).is_ok_and(|x| {
        close(x[0].re, T::from_f64(2.0), 8.0)
            && close(x[0].eps, T::from_f64(-4.0), 8.0)
            && close(x[1].re, T::from_f64(-1.0), 8.0)
            && close(x[1].eps, T::from_f64(2.0), 8.0)
    }));
}

/// Sehr kleine und sehr große Einträge: Ihr Quadrat läuft unter bzw. über,
/// die Matrizen sind trotzdem regulär und die Normen darstellbar.
fn extremes<T: Real>(tiny: f64, huge: f64) {
    let (tiny, huge) = (T::from_f64(tiny), T::from_f64(huge));
    for value in [tiny, huge] {
        assert!(Matrix([[value]]).solve(Vector([value])).is_ok_and(|x| close(x[0], T::ONE, 8.0)));
        assert!(Matrix([[value, T::ZERO], [T::ZERO, value]])
            .inverse()
            .is_ok_and(|inverse| close(inverse[(0, 0)] * value, T::ONE, 8.0)));
        let z = Complex::new(T::ZERO, value);
        assert!(Matrix([[z]]).lu().is_ok());
    }
    assert!(matches!(Matrix([[T::ZERO]]).lu(), Err(LinalgError::Singular)));
    assert!(matches!(Matrix([[T::NAN]]).lu(), Err(LinalgError::Singular)));
}

#[test]
fn extremes_f32() {
    extremes::<f32>(1.0e-30, 1.0e30);
}

#[test]
fn extremes_f64() {
    extremes::<f64>(1.0e-200, 1.0e200);
}

#[test]
fn lu_f32() {
    lu::<f32>();
}

#[test]
fn lu_f64() {
    lu::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod with_roots {
    use super::*;
    use crate::scalar::RealElementary;

    fn close_cm<T: Real, const R: usize, const C: usize>(
        a: Matrix<Complex<T>, R, C>,
        b: Matrix<Complex<T>, R, C>,
        ulps: f64,
    ) -> bool {
        close_m(a.map(|z| z.re), b.map(|z| z.re), ulps) && close_m(a.map(|z| z.im), b.map(|z| z.im), ulps)
    }

    fn cholesky<T: RealElementary>() {
        let id2 = Matrix::<T, 2, 2>::identity();
        let a = m::<T, 2, 2>([[4.0, 2.0], [2.0, 3.0]]);
        let root2 = T::from_f64(2.0).sqrt();
        let expected = Matrix([[T::from_f64(2.0), T::ZERO], [T::ONE, root2]]);

        assert!(a.cholesky().is_ok_and(|ch| {
            let l = ch.l();
            close_m(l, expected, 8.0)
                && close_m(l * l.adjoint(), a, 8.0)
                // 4x + 2y = 2, 2x + 3y = 1  →  x = 0.5, y = 0
                && close_v(ch.solve(v([2.0, 1.0])), v([0.5, 0.0]), 8.0)
                && close(ch.det(), T::from_f64(8.0), 8.0)
                && close_m(ch.inverse() * a, id2, 8.0)
        }));

        // Nur das untere Dreieck zählt.
        let upper_garbage = m::<T, 2, 2>([[4.0, 99.0], [2.0, 3.0]]);
        assert!(upper_garbage.cholesky().is_ok_and(|ch| close_m(ch.l(), expected, 8.0)));

        // Nicht positiv definit.
        assert!(matches!(
            m::<T, 2, 2>([[1.0, 2.0], [2.0, 1.0]]).cholesky(),
            Err(LinalgError::NotPositiveDefinite)
        ));
        assert!(matches!(
            Matrix::<T, 2, 2>::ZERO.cholesky(),
            Err(LinalgError::NotPositiveDefinite)
        ));

        // Komplex hermitesch: [[2, i], [−i, 2]] hat die Eigenwerte 1 und 3.
        let i = Complex::<T>::I;
        let two = Complex::from_real(T::from_f64(2.0));
        let h = Matrix([[two, i], [-i, two]]);
        assert!(h.cholesky().is_ok_and(|ch| {
            let l = ch.l();
            let x = ch.solve(Vector([two, Complex::ZERO]));
            let residual = h * x;
            close_cm(l * l.adjoint(), h, 8.0)
                && close(residual[0].re, T::from_f64(2.0), 8.0)
                && close(residual[0].im, T::ZERO, 8.0)
                && close(residual[1].re, T::ZERO, 8.0)
                && close(residual[1].im, T::ZERO, 8.0)
        }));

        // Ableitung: L = √p, dL/dp = 1/(2√p); bei p = 4: (2, 0.25).
        let dual = Matrix([[Dual::variable(T::from_f64(4.0))]]);
        assert!(dual.cholesky().is_ok_and(|ch| {
            let l = ch.l()[(0, 0)];
            close(l.re, T::from_f64(2.0), 8.0) && close(l.eps, T::from_f64(0.25), 8.0)
        }));
    }

    fn qr<T: RealElementary>() {
        // Überbestimmt: 3 Gleichungen, 2 Unbekannte.
        let a = m::<T, 3, 2>([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]);
        let dec = a.qr();
        let (q, r) = (dec.q(), dec.r());
        assert!(close_m(q * r, a, 8.0));
        assert!(close_m(q.adjoint() * q, Matrix::identity(), 8.0));
        assert_eq!((r[(1, 0)], r[(2, 0)], r[(2, 1)]), (T::ZERO, T::ZERO, T::ZERO));

        // Exakt lösbar: b = A·(1, 1).
        assert!(dec.solve(v([1.0, 1.0, 2.0])).is_ok_and(|x| close_v(x, v([1.0, 1.0]), 16.0)));
        // Nicht lösbar: Normalgleichungen [[2, 1], [1, 2]]·x = (1, 1)
        // ergeben x = (1/3, 1/3).
        let third = T::ONE / T::from_f64(3.0);
        assert!(dec
            .solve(v([1.0, 1.0, 0.0]))
            .is_ok_and(|x| close_v(x, Vector([third, third]), 16.0)));

        // Quadratisch, Diagonalelement null.
        let swap = m::<T, 2, 2>([[0.0, 1.0], [1.0, 0.0]]);
        let dec = swap.qr();
        assert!(close_m(dec.q() * dec.r(), swap, 8.0));
        assert!(dec.solve(v([2.0, 3.0])).is_ok_and(|x| close_v(x, v([3.0, 2.0]), 8.0)));

        // Rangabfall und zu wenige Zeilen.
        let deficient = m::<T, 3, 2>([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0]]);
        assert_eq!(deficient.qr().solve(v([1.0, 1.0, 1.0])), Err(LinalgError::Singular));
        let wide = m::<T, 1, 2>([[1.0, 2.0]]);
        assert_eq!(wide.qr().solve(v([1.0])), Err(LinalgError::DimensionMismatch));
        assert!(close_m(wide.qr().q() * wide.qr().r(), wide, 8.0));

        // Komplex.
        let i = Complex::<T>::I;
        let one = Complex::<T>::ONE;
        let z = Matrix([[i, one], [one, i]]);
        let dec = z.qr();
        assert!(close_cm(dec.q() * dec.r(), z, 8.0));
        assert!(close_cm(dec.q().adjoint() * dec.q(), Matrix::identity(), 8.0));
        assert!(dec.solve(Vector([one, Complex::ZERO])).is_ok_and(|x| {
            close(x[0].re, T::ZERO, 8.0)
                && close(x[0].im, T::from_f64(-0.5), 8.0)
                && close(x[1].re, T::from_f64(0.5), 8.0)
                && close(x[1].im, T::ZERO, 8.0)
        }));

        // Ableitung durch die Ausgleichsrechnung: p·x = 1 → x = 1/p,
        // dx/dp = −1/p²; bei p = 2: (0.5, −0.25).
        let dual = Matrix([[Dual::variable(T::from_f64(2.0))]]);
        assert!(dual.qr().solve(Vector([Dual::constant(T::ONE)])).is_ok_and(|x| {
            close(x[0].re, T::from_f64(0.5), 8.0) && close(x[0].eps, T::from_f64(-0.25), 8.0)
        }));
    }

    #[test]
    fn cholesky_f32() {
        cholesky::<f32>();
    }

    #[test]
    fn cholesky_f64() {
        cholesky::<f64>();
    }

    #[test]
    fn qr_f32() {
        qr::<f32>();
    }

    #[test]
    fn qr_f64() {
        qr::<f64>();
    }
}
