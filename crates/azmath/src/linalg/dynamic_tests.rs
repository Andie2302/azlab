//! Tests für Sichten, `LinearOperator` und CSR-Matrizen.
//!
//! Geprüft werden Struktur, Fehlerfälle und die Übereinstimmung mit den
//! Matrizen fester Größe; es gibt keine vorab ausgerechneten Erwartungswerte.

use super::{CsrMatrix, LinalgError, LinearOperator, Matrix, MatrixView, MatrixViewMut, Vector};
use crate::scalar::{Algebra, Real};
use crate::testutil::close;

fn check<T: Real>() {
    let f = T::from_f64;
    let data = [f(4.0), f(1.0), f(-0.5), f(-1.0), f(5.0), f(0.8), f(0.3), f(-0.7), f(3.0)];
    let fixed = Matrix([
        [data[0], data[1], data[2]],
        [data[3], data[4], data[5]],
        [data[6], data[7], data[8]],
    ]);
    let b = [f(1.0), f(-2.0), f(0.5)];

    // Sicht: Größe und Lage der Elemente im Puffer.
    assert!(MatrixView::new(&data, 3, 3).is_ok_and(|view| {
        view.rows() == 3
            && view.cols() == 3
            && view[(1, 2)] == data[5]
            && view.get(2, 0) == Some(&data[6])
            && view.get(3, 0).is_none()
            && view.get(0, 3).is_none()
            && view.data().len() == 9
    }));
    // Dieselben Daten als 2 × 4: Die ersten acht Elemente, anders aufgeteilt.
    let eight = data.get(..8).unwrap_or(&[]);
    assert!(MatrixView::new(eight, 2, 4).is_ok_and(|view| view[(1, 0)] == data[4]));
    assert!(matches!(MatrixView::new(&data, 2, 3), Err(LinalgError::DimensionMismatch)));

    let mut buffer = data;
    assert!(matches!(
        MatrixViewMut::new(&mut buffer, 4, 4),
        Err(LinalgError::DimensionMismatch)
    ));
    assert!(MatrixViewMut::new(&mut buffer, 3, 3).is_ok_and(|mut view| {
        view[(0, 0)] = T::ONE;
        if let Some(slot) = view.get_mut(2, 2) {
            *slot = T::ZERO;
        }
        view.rows() == 3
            && view.cols() == 3
            && view[(0, 0)] == T::ONE
            && view.get(2, 2) == Some(&T::ZERO)
            && view.get(3, 3).is_none()
            && view.as_view()[(0, 1)] == data[1]
    }));
    assert_eq!(buffer[0], T::ONE);
    assert_eq!(buffer[8], T::ZERO);

    // LinearOperator: alle Speicherformen wirken wie Matrix · Vektor.
    let expected = fixed * Vector(b);
    let agrees = |operator: &dyn LinearOperator<T>| {
        let mut y = [T::ZERO; 3];
        operator.rows() == 3
            && operator.cols() == 3
            && operator.apply(&b, &mut y).is_ok()
            && y.iter().zip(expected.0).all(|(&got, want)| close(got, want, 8.0))
    };
    assert!(agrees(&fixed));
    assert!(MatrixView::new(&data, 3, 3).is_ok_and(|view| agrees(&view)));
    let mut copy = data;
    assert!(MatrixViewMut::new(&mut copy, 3, 3).is_ok_and(|view| agrees(&view)));

    let col_idx = [0, 1, 2, 0, 1, 2, 0, 1, 2];
    let row_ptr = [0, 3, 6, 9];
    assert!(CsrMatrix::new(3, 3, &row_ptr, &col_idx, &data).is_ok_and(|csr| {
        agrees(&csr)
            && csr.nnz() == 9
            && csr.get(1, 2) == Some(data[5])
            && csr.get(3, 0).is_none()
            && csr.get(0, 3).is_none()
    }));

    // Wirklich dünnbesetzt: nur die Diagonale gespeichert.
    let diagonal = [data[0], data[4], data[8]];
    assert!(CsrMatrix::new(3, 3, &[0, 1, 2, 3], &[0, 1, 2], &diagonal).is_ok_and(|csr| {
        let mut y = [T::ZERO; 3];
        csr.apply(&b, &mut y).is_ok()
            && y.iter().zip(diagonal.iter().zip(b)).all(|(&got, (&d, x))| got == d * x)
            && csr.get(0, 1) == Some(T::ZERO)
    }));

    // Falsche Längen.
    let mut y = [T::ZERO; 2];
    assert_eq!(fixed.apply(&b, &mut y), Err(LinalgError::DimensionMismatch));
    assert!(MatrixView::new(&data, 3, 3)
        .is_ok_and(|view| view.apply(&b, &mut y) == Err(LinalgError::DimensionMismatch)));
    assert!(CsrMatrix::new(3, 3, &row_ptr, &col_idx, &data)
        .is_ok_and(|csr| csr.apply(&b, &mut y) == Err(LinalgError::DimensionMismatch)));

    // Fehlerhafte CSR-Strukturen.
    let broken = |row_ptr: &[usize], col_idx: &[usize]| {
        matches!(
            CsrMatrix::new(3, 3, row_ptr, col_idx, &diagonal),
            Err(LinalgError::DimensionMismatch)
        )
    };
    assert!(broken(&[0, 1, 2], &[0, 1, 2]));
    assert!(broken(&[1, 1, 2, 3], &[0, 1, 2]));
    assert!(broken(&[0, 2, 1, 3], &[0, 1, 2]));
    assert!(broken(&[0, 1, 2, 4], &[0, 1, 2]));
    assert!(broken(&[0, 1, 2, 3], &[0, 1]));
    assert!(broken(&[0, 1, 2, 3], &[0, 1, 3]));

    // LU auf einer Sicht liefert dieselbe Lösung wie LU fester Größe.
    let mut factors = data;
    let mut perm = [0; 3];
    let mut x = [T::ZERO; 3];
    assert!(MatrixViewMut::new(&mut factors, 3, 3).is_ok_and(|mut view| {
        view.lu_in_place(&mut perm).is_ok() && view.as_view().lu_solve(&perm, &b, &mut x).is_ok()
    }));
    assert!(fixed
        .solve(Vector(b))
        .is_ok_and(|exact| x.iter().zip(exact.0).all(|(&got, want)| close(got, want, 64.0))));

    // Fehlerfälle der Zerlegung.
    let mut rectangular = [T::ONE; 6];
    assert!(MatrixViewMut::new(&mut rectangular, 2, 3).is_ok_and(|mut view| {
        view.lu_in_place(&mut perm) == Err(LinalgError::DimensionMismatch)
    }));
    let mut singular = [T::ONE; 4];
    let mut pair = [0; 2];
    assert!(MatrixViewMut::new(&mut singular, 2, 2)
        .is_ok_and(|mut view| view.lu_in_place(&mut pair) == Err(LinalgError::Singular)));
    let mut wrong_perm = [0; 2];
    let mut again = data;
    assert!(MatrixViewMut::new(&mut again, 3, 3).is_ok_and(|mut view| {
        view.lu_in_place(&mut wrong_perm) == Err(LinalgError::DimensionMismatch)
    }));
    assert_eq!(<Matrix<T, 3, 3> as Algebra>::ONE, Matrix::identity());
}

#[test]
fn dynamic_f32() {
    check::<f32>();
}

#[test]
fn dynamic_f64() {
    check::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod with_roots {
    use super::*;
    use crate::scalar::RealElementary;

    fn cholesky<T: RealElementary>() {
        let f = T::from_f64;
        // Symmetrisch positiv definit durch Konstruktion: MᵀM + I.
        let m = Matrix([[f(1.0), f(0.5)], [f(0.2), f(1.5)]]);
        let spd = m.transpose() * m + Matrix::identity();
        let mut buffer = [spd[(0, 0)], spd[(0, 1)], spd[(1, 0)], spd[(1, 1)]];
        assert!(MatrixViewMut::new(&mut buffer, 2, 2)
            .is_ok_and(|mut view| view.cholesky_in_place().is_ok()));
        // Dieselben Faktoren wie bei der Matrix fester Größe.
        assert!(spd.cholesky().is_ok_and(|fixed| {
            let l = fixed.l();
            [l[(0, 0)], l[(0, 1)], l[(1, 0)], l[(1, 1)]]
                .into_iter()
                .zip(buffer)
                .all(|(want, got)| close(got, want, 8.0))
        }));

        let mut rectangular = [T::ONE; 6];
        assert!(MatrixViewMut::new(&mut rectangular, 2, 3).is_ok_and(|mut view| {
            view.cholesky_in_place() == Err(LinalgError::DimensionMismatch)
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
}

#[cfg(feature = "alloc")]
mod owned {
    use super::*;
    use crate::linalg::DMatrix;
    use alloc::vec;

    fn check<T: Real>() {
        let f = T::from_f64;
        let a = Matrix([[f(4.0), f(1.0), f(-0.5)], [f(-1.0), f(5.0), f(0.8)], [f(0.3), f(-0.7), f(3.0)]]);
        let c = Matrix([[f(1.0), f(0.5)], [f(0.2), f(1.5)], [f(-0.6), f(0.1)]]);
        let b = [f(1.0), f(-2.0), f(0.5)];

        let da = DMatrix::from_fn(3, 3, |i, j| a[(i, j)]);
        let dc = DMatrix::from_fn(3, 2, |i, j| c[(i, j)]);
        let ok = match (da, dc) {
            (Ok(da), Ok(dc)) => {
                let fixed_product = a * c;
                let product_ok = da.mul(&dc).is_ok_and(|p| {
                    p.rows() == 3
                        && p.cols() == 2
                        && (0..3).all(|i| (0..2).all(|j| p[(i, j)] == fixed_product[(i, j)]))
                }) && dc.mul(&da).is_err();

                let vector_ok = da.mul_vec(&b).is_ok_and(|y| y == (a * Vector(b)).0)
                    && da.mul_vec(&b[..2]).is_err();

                let transpose_ok = {
                    let t = dc.transpose();
                    t.rows() == 2
                        && t.cols() == 3
                        && (0..3).all(|i| (0..2).all(|j| t[(j, i)] == dc[(i, j)]))
                        && t.transpose() == dc
                };

                let solve_ok = match (da.solve(&b), a.solve(Vector(b))) {
                    (Ok(x), Ok(exact)) => {
                        x.iter().zip(exact.0).all(|(&got, want)| close(got, want, 64.0))
                    }
                    _ => false,
                } && dc.solve(&b).is_err()
                    && da.solve(&b[..2]).is_err();

                // LinearOperator und Sichten.
                let mut y = [T::ZERO; 3];
                let operator_ok = da.apply(&b, &mut y).is_ok()
                    && y == (a * Vector(b)).0
                    && LinearOperator::rows(&da) == 3
                    && LinearOperator::cols(&dc) == 2
                    && da.view()[(1, 2)] == a[(1, 2)]
                    && da.get(1, 2) == Some(&a[(1, 2)])
                    && da.get(3, 0).is_none();
                let mut copy = da.clone();
                copy.view_mut()[(0, 0)] = T::ZERO;
                copy[(1, 1)] = T::ONE;
                let mutation_ok = copy[(0, 0)] == T::ZERO && copy[(1, 1)] == T::ONE && copy != da;

                product_ok && vector_ok && transpose_ok && solve_ok && operator_ok && mutation_ok
            }
            _ => false,
        };
        assert!(ok);

        // Erzeugen.
        assert!(DMatrix::<T>::identity(3).is_ok_and(|id| {
            (0..3).all(|i| (0..3).all(|j| id[(i, j)] == Matrix::<T, 3, 3>::identity()[(i, j)]))
        }));
        assert!(DMatrix::<T>::zeros(2, 3).is_ok_and(|z| z.rows() == 2 && z.cols() == 3));
        assert!(DMatrix::<T>::zeros(usize::MAX, 2).is_err());
        assert!(DMatrix::from_vec(vec![T::ONE; 6], 2, 3).is_ok());
        assert!(DMatrix::from_vec(vec![T::ONE; 6], 2, 2).is_err());
        // Singulär.
        assert!(DMatrix::from_vec(vec![T::ONE; 4], 2, 2)
            .is_ok_and(|m| m.solve(&[T::ONE, T::ONE]) == Err(LinalgError::Singular)));
    }

    #[test]
    fn owned_f32() {
        check::<f32>();
    }

    #[test]
    fn owned_f64() {
        check::<f64>();
    }
}
