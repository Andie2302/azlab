//! Tests für Eigenwerte und Singulärwerte.

#![cfg(any(feature = "std", feature = "libm"))]

use crate::linalg::{LinalgError, Matrix, Vector};
use crate::number::Complex;
use crate::scalar::{Algebra, RealElementary};
use crate::testutil::close;

const ULPS: f64 = 2048.0;

fn m<T: RealElementary, const R: usize, const C: usize>(rows: [[f64; C]; R]) -> Matrix<T, R, C> {
    Matrix(rows.map(|row| row.map(T::from_f64)))
}

fn close_m<T: RealElementary, const R: usize, const C: usize>(
    a: Matrix<T, R, C>,
    b: Matrix<T, R, C>,
) -> bool {
    (0..R).all(|i| (0..C).all(|j| close(a[(i, j)], b[(i, j)], ULPS)))
}

/// Jeder erwartete Wert kommt unter den berechneten vor.
fn contains_all<T: RealElementary>(found: &[Complex<T>], expected: &[(f64, f64)]) -> bool {
    let limit = T::EPSILON.sqrt();
    expected.iter().all(|&(re, im)| {
        let target = Complex::new(T::from_f64(re), T::from_f64(im));
        found.iter().any(|&value| (value - target).norm() <= limit)
    })
}

fn symmetric<T: RealElementary>() {
    let f = T::from_f64;
    // [[2, 1], [1, 2]]: Eigenwerte 1 und 3, Eigenvektoren (1, ∓1)/√2.
    let eigen = m::<T, 2, 2>([[2.0, 1.0], [1.0, 2.0]]).symmetric_eigen();
    assert!(eigen.is_ok_and(|e| {
        let half = f(0.5).sqrt();
        close(e.values[0], f(1.0), ULPS)
            && close(e.values[1], f(3.0), ULPS)
            && close(e.vectors[(0, 0)].abs(), half, ULPS)
            && close(e.vectors[(1, 1)].abs(), half, ULPS)
    }));

    let a = m::<T, 4, 4>([
        [4.0, 1.0, -2.0, 0.5],
        [1.0, 3.0, 0.0, -1.0],
        [-2.0, 0.0, 5.0, 2.0],
        [0.5, -1.0, 2.0, 1.0],
    ]);
    assert!(a.symmetric_eigen().is_ok_and(|e| {
        let v = e.vectors;
        // Aufsteigend, orthonormal, A·V = V·D und A = V·D·Vᵀ, Spur erhalten.
        e.values.0.windows(2).all(|pair| pair[0] <= pair[1])
            && close_m(v.transpose() * v, Matrix::identity())
            && close_m(a * v, v * Matrix::from_diagonal(e.values))
            && close_m(v * Matrix::from_diagonal(e.values) * v.transpose(), a)
            && close(e.values.sum(), a.trace(), ULPS)
    }));
    // Schon diagonal, null, 1×1.
    assert!(m::<T, 3, 3>([[3.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 2.0]])
        .symmetric_eigen()
        .is_ok_and(|e| e.values == Vector([f(-1.0), f(2.0), f(3.0)])));
    assert!(Matrix::<T, 3, 3>::ZERO.symmetric_eigen().is_ok_and(|e| e.values == Vector::ZERO));
    assert!(m::<T, 1, 1>([[7.0]]).symmetric_eigen().is_ok_and(|e| e.values[0] == f(7.0)));
    assert_eq!(
        Matrix([[T::NAN, T::ONE], [T::ONE, T::ONE]]).symmetric_eigen().err(),
        Some(LinalgError::NotConverged)
    );
}

fn general<T: RealElementary>() {
    let f = T::from_f64;
    // Drehung um 90°: Eigenwerte ±i.
    let rotation = m::<T, 2, 2>([[0.0, -1.0], [1.0, 0.0]]);
    assert!(rotation.eigenvalues().is_ok_and(|v| contains_all(&v, &[(0.0, 1.0), (0.0, -1.0)])));
    // Begleitmatrix von (x − 1)(x − 2)(x − 3) = x³ − 6x² + 11x − 6.
    let companion = m::<T, 3, 3>([[0.0, 0.0, 6.0], [1.0, 0.0, -11.0], [0.0, 1.0, 6.0]]);
    assert!(companion
        .eigenvalues()
        .is_ok_and(|v| contains_all(&v, &[(1.0, 0.0), (2.0, 0.0), (3.0, 0.0)])));
    // Dreiecksmatrix: Eigenwerte auf der Diagonale.
    let triangular = m::<T, 3, 3>([[2.0, 5.0, -1.0], [0.0, -3.0, 4.0], [0.0, 0.0, 0.5]]);
    assert!(triangular
        .eigenvalues()
        .is_ok_and(|v| contains_all(&v, &[(2.0, 0.0), (-3.0, 0.0), (0.5, 0.0)])));
    // Allgemein: Summe = Spur, Produkt = Determinante; reelle Matrix →
    // Imaginärteile heben sich auf.
    let a = m::<T, 4, 4>([
        [1.0, 2.0, 0.0, -1.0],
        [-3.0, 0.5, 1.0, 2.0],
        [0.0, -1.0, 2.0, 1.5],
        [2.0, 0.0, -2.0, 1.0],
    ]);
    let limit = T::EPSILON.sqrt();
    assert!(a.eigenvalues().is_ok_and(|v| {
        let sum = v.iter().fold(Complex::ZERO, |acc, &value| acc + value);
        let product = v.iter().fold(Complex::ONE, |acc, &value| acc * value);
        (sum.re - a.trace()).abs() <= limit
            && sum.im.abs() <= limit
            && (product.re - a.det()).abs() <= limit * f(10.0)
            && product.im.abs() <= limit * f(10.0)
    }));
    // Symmetrisch: dieselben Werte wie das Jacobi-Verfahren.
    let s = m::<T, 3, 3>([[2.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 2.0]]);
    assert!(match (s.eigenvalues(), s.symmetric_eigen()) {
        (Ok(general), Ok(jacobi)) => jacobi.values.iter().all(|&value| {
            general.iter().any(|z| (z.re - value).abs() <= limit && z.im.abs() <= limit)
        }),
        _ => false,
    });
    // Komplexe Matrix.
    let i = Complex::<T>::I;
    let two = Complex::from_real(f(2.0));
    let z = Matrix([[i, Complex::ONE], [Complex::ZERO, two]]);
    assert!(z.eigenvalues().is_ok_and(|v| contains_all(&v, &[(0.0, 1.0), (2.0, 0.0)])));
    let full = Matrix([[i, Complex::ONE], [Complex::ONE, i]]);
    // Eigenwerte i ± 1.
    assert!(full.eigenvalues().is_ok_and(|v| contains_all(&v, &[(1.0, 1.0), (-1.0, 1.0)])));
    // Sonderfälle.
    assert!(m::<T, 1, 1>([[4.0]]).eigenvalues().is_ok_and(|v| v[0] == Complex::from_real(f(4.0))));
    assert!(Matrix::<T, 0, 0>::ZERO.eigenvalues().is_ok());
    assert!(Matrix::<T, 3, 3>::ZERO.eigenvalues().is_ok_and(|v| contains_all(&v, &[(0.0, 0.0)])));
    assert_eq!(
        Matrix([[T::NAN, T::ONE], [T::ONE, T::ONE]]).eigenvalues().err(),
        Some(LinalgError::NotConverged)
    );
}

fn singular_values<T: RealElementary>() {
    let f = T::from_f64;
    // Diagonal: Singulärwerte sind die Beträge, absteigend.
    assert!(m::<T, 2, 2>([[3.0, 0.0], [0.0, -4.0]])
        .svd()
        .is_ok_and(|d| close(d.s[0], f(4.0), ULPS) && close(d.s[1], f(3.0), ULPS)));

    let a = m::<T, 4, 3>([
        [1.0, 2.0, 0.0],
        [-3.0, 0.5, 1.0],
        [0.0, -1.0, 2.0],
        [2.0, 0.0, -2.0],
    ]);
    assert!(a.svd().is_ok_and(|d| {
        d.s.0.windows(2).all(|pair| pair[0] >= pair[1])
            && d.s.iter().all(|&value| value >= T::ZERO)
            && close_m(d.u.transpose() * d.u, Matrix::identity())
            && close_m(d.v.transpose() * d.v, Matrix::identity())
            && close_m(d.u * Matrix::from_diagonal(d.s) * d.v.transpose(), a)
            && d.rank() == 3
            && d.condition_number() >= T::ONE
            // Voller Spaltenrang: A⁺·A = I.
            && close_m(d.pseudo_inverse() * a, Matrix::identity())
    }));
    // Singulärwerte sind die Wurzeln der Eigenwerte von AᵀA.
    assert!(match (a.svd(), (a.transpose() * a).symmetric_eigen()) {
        (Ok(d), Ok(e)) => (0..3).all(|k| close(d.s[k] * d.s[k], e.values[2 - k], ULPS)),
        _ => false,
    });
    // Ausgleichsrechnung: dieselbe Lösung wie über QR.
    let b = Vector([f(1.0), f(-2.0), f(0.5), f(3.0)]);
    assert!(match (a.svd(), a.qr().solve(b)) {
        (Ok(d), Ok(x)) => d.solve(b).iter().zip(x.iter()).all(|(&p, &q)| close(p, q, ULPS)),
        _ => false,
    });

    // Rangabfall: zweite Spalte ist das Doppelte der ersten.
    let deficient = m::<T, 3, 2>([[1.0, 2.0], [2.0, 4.0], [-1.0, -2.0]]);
    assert!(deficient.svd().is_ok_and(|d| {
        let pinv = d.pseudo_inverse();
        d.rank() == 1
            && d.s[1] <= T::EPSILON.sqrt()
            // Eigenschaft der Pseudoinversen: A·A⁺·A = A.
            && close_m(deficient * pinv * deficient, deficient)
            && d.condition_number() > f(1.0e3)
    }));
    assert!(Matrix::<T, 2, 2>::ZERO.svd().is_ok_and(|d| {
        d.rank() == 0 && d.s == Vector::ZERO && d.condition_number() == T::INFINITY
    }));
    assert_eq!(m::<T, 1, 2>([[1.0, 2.0]]).svd().err(), Some(LinalgError::DimensionMismatch));
    assert_eq!(
        Matrix([[T::NAN, T::ONE], [T::ONE, T::ONE]]).svd().err(),
        Some(LinalgError::NotConverged)
    );
}

/// Zerlegung über viele Größenordnungen: Rekonstruktion und
/// Orthogonalität müssen unabhängig vom Maßstab stimmen.
fn magnitudes<T: RealElementary>(scales: [f64; 4]) {
    let unit = m::<T, 2, 2>([[1.0, 1.0], [0.0, 1.0]]);
    let reference = unit.svd();
    for scale in scales.map(T::from_f64) {
        let a = unit * scale;
        assert!(match (a.svd(), reference) {
            (Ok(d), Ok(r)) => {
                close_m(d.u.transpose() * d.u, Matrix::identity())
                    && close_m(d.v.transpose() * d.v, Matrix::identity())
                    // Rekonstruktion, gemessen am Maßstab.
                    && close_m((d.u * Matrix::from_diagonal(d.s) * d.v.transpose()) / scale, unit)
                    && close(d.s[0] / scale, r.s[0], ULPS)
                    && close(d.s[1] / scale, r.s[1], ULPS)
                    && d.rank() == 2
                    && close(d.condition_number(), r.condition_number(), ULPS)
                    && close_m(d.pseudo_inverse() * a, Matrix::identity())
            }
            _ => false,
        });
    }
}

/// Verschiedene Größenordnungen innerhalb derselben Matrix: Die Quadrate der
/// kleinen Spalten laufen unter, während die große Spalte den gemeinsamen
/// Maßstab bestimmt.
fn mixed_magnitudes<T: RealElementary>(small: f64) {
    let f = T::from_f64;
    let small = f(small);
    assert!(small > T::ZERO && small * small == T::ZERO);
    // Block diag(1, small·[[1, 1], [0, 1]]): Die Singulärwerte des kleinen
    // Blocks sind small·(√5 ± 1)/2.
    let golden = (f(5.0).sqrt() + T::ONE) * f(0.5);
    let expected = [T::ONE, small * golden, small * (golden - T::ONE)];
    for (big, shift) in [(0, 1), (1, 2), (2, 0)] {
        // Dieselbe Matrix mit vertauschter Reihenfolge der Spalten und
        // Zeilen: Die große Spalte steht vorn, in der Mitte oder hinten.
        let mut a = Matrix::<T, 3, 3>::ZERO;
        let (first, second) = ((big + 1) % 3, (big + 2) % 3);
        a[(shift, big)] = T::ONE;
        a[((shift + 1) % 3, first)] = small;
        a[((shift + 1) % 3, second)] = small;
        a[((shift + 2) % 3, second)] = small;
        assert!(a.svd().is_ok_and(|d| {
            let column = |matrix: Matrix<T, 3, 3>, j: usize| Vector::<T, 3>(core::array::from_fn(|i| matrix[(i, j)]));
            close_m(d.u.transpose() * d.u, Matrix::identity())
                && close_m(d.v.transpose() * d.v, Matrix::identity())
                && (0..3).all(|j| close(d.s[j], expected[j], ULPS))
                // A·v = s·u für jeden Singulärwert, gemessen an seiner
                // eigenen Größe.
                && (0..3).all(|j| {
                    let (image, target) = (a * column(d.v, j), column(d.u, j) * d.s[j]);
                    (image - target).norm() <= d.s[j] * T::EPSILON * f(ULPS)
                })
        }));
    }
}

#[test]
fn magnitudes_f32() {
    magnitudes::<f32>([1.0e15, 1.0e-15, 1.0e25, 1.0e-25]);
    mixed_magnitudes::<f32>(1.0e-25);
}

#[test]
fn magnitudes_f64() {
    magnitudes::<f64>([1.0e100, 1.0e-100, 1.0e200, 1.0e-200]);
    mixed_magnitudes::<f64>(1.0e-200);
}

#[test]
fn spectral_f32() {
    symmetric::<f32>();
    general::<f32>();
    singular_values::<f32>();
}

#[test]
fn spectral_f64() {
    symmetric::<f64>();
    general::<f64>();
    singular_values::<f64>();
}
