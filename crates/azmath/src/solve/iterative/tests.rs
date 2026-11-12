//! Tests der iterativen Löser: Vergleich mit der LU-Lösung und Restfehler.

#![cfg(any(feature = "std", feature = "libm"))]

use super::{bicgstab, conjugate_gradient, gmres, gmres_work_len, preconditioned_conjugate_gradient};
use crate::linalg::{CsrMatrix, LinearOperator, Matrix, Vector};
use crate::scalar::RealElementary;
use crate::solve::{SolveError, Tolerance};

/// Abstand zweier Vektoren, gemessen an der Wurzel der Maschinengenauigkeit.
fn near<T: RealElementary>(a: &[T], b: &[T]) -> bool {
    let limit = T::EPSILON.sqrt() * T::from_f64(64.0);
    a.len() == b.len() && a.iter().zip(b).all(|(&x, &y)| (x - y).abs() <= limit * (T::ONE + y.abs()))
}

fn check<T: RealElementary>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default()
        .with_rel_x(T::EPSILON.sqrt())
        .with_max_iterations(200);
    let b = [f(1.0), f(-2.0), f(0.5)];

    // Symmetrisch positiv definit durch Konstruktion: MᵀM + I.
    let m = Matrix([
        [f(1.0), f(0.5), f(-0.3)],
        [f(0.2), f(1.5), f(0.4)],
        [f(-0.6), f(0.1), f(2.0)],
    ]);
    let spd = m.transpose() * m + Matrix::identity();
    let exact = spd.solve(Vector(b));
    for solver in [conjugate_gradient::<T, Matrix<T, 3, 3>>, bicgstab::<T, Matrix<T, 3, 3>>] {
        let mut x = [T::ZERO; 3];
        let mut work = [T::ZERO; 18];
        assert!(solver(&spd, &b, &mut x, &mut work, &tol).is_ok());
        assert!(exact.is_ok_and(|e| near(&x, &e.0)));
        // Rückeinsetzen.
        let mut image = [T::ZERO; 3];
        assert!(spd.apply(&x, &mut image).is_ok() && near(&image, &b));
    }

    // Unsymmetrisch, diagonal dominant: nur BiCGSTAB.
    let general = Matrix([
        [f(4.0), f(1.0), f(-0.5)],
        [f(-1.0), f(5.0), f(0.8)],
        [f(0.3), f(-0.7), f(3.0)],
    ]);
    let mut x = [T::ZERO; 3];
    let mut work = [T::ZERO; 18];
    assert!(bicgstab(&general, &b, &mut x, &mut work, &tol).is_ok());
    assert!(general.solve(Vector(b)).is_ok_and(|e| near(&x, &e.0)));

    // Dieselbe Matrix dünnbesetzt gespeichert liefert dieselbe Lösung.
    let values = [
        f(4.0), f(1.0), f(-0.5), f(-1.0), f(5.0), f(0.8), f(0.3), f(-0.7), f(3.0),
    ];
    let col_idx = [0, 1, 2, 0, 1, 2, 0, 1, 2];
    let row_ptr = [0, 3, 6, 9];
    let sparse = CsrMatrix::new(3, 3, &row_ptr, &col_idx, &values);
    assert!(sparse.is_ok_and(|csr| {
        let mut y = [T::ZERO; 3];
        let mut work = [T::ZERO; 18];
        bicgstab(&csr, &b, &mut y, &mut work, &tol).is_ok() && near(&y, &x)
    }));

    // Die Lösung als Startwert: sofort fertig.
    let mut work = [T::ZERO; 18];
    let mut again = x;
    assert!(bicgstab(&general, &b, &mut again, &mut work, &tol).is_ok_and(|s| s.iterations == 0));

    // Fehlerfälle.
    let mut x = [T::ZERO; 3];
    let mut small = [T::ZERO; 5];
    assert_eq!(
        conjugate_gradient(&spd, &b, &mut x, &mut small, &tol),
        Err(SolveError::InvalidInput)
    );
    assert_eq!(bicgstab(&spd, &b, &mut x, &mut small, &tol), Err(SolveError::InvalidInput));
    let mut work = [T::ZERO; 18];
    let mut short = [T::ZERO; 2];
    assert_eq!(
        conjugate_gradient(&spd, &b, &mut short, &mut work, &tol),
        Err(SolveError::InvalidInput)
    );
    // Nicht positiv definit: CG meldet es.
    let indefinite = Matrix([[f(1.0), f(0.0), f(0.0)], [f(0.0), f(-1.0), f(0.0)], [f(0.0), f(0.0), f(1.0)]]);
    assert!(matches!(
        conjugate_gradient(&indefinite, &[f(0.0), f(1.0), f(0.0)], &mut x, &mut work, &tol),
        Err(SolveError::Singular(_))
    ));
    // Zu wenige Iterationen.
    let none = tol.with_max_iterations(0);
    let mut x = [T::ZERO; 3];
    assert!(matches!(
        conjugate_gradient(&spd, &b, &mut x, &mut work, &none),
        Err(SolveError::NotConverged(s)) if s.iterations == 0
    ));
}

fn check_gmres<T: RealElementary>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default()
        .with_rel_x(T::EPSILON.sqrt())
        .with_max_iterations(200);
    let b = [f(1.0), f(-2.0), f(0.5), f(0.25)];
    // Unsymmetrisch, diagonal dominant.
    let a = Matrix([
        [f(4.0), f(1.0), f(-0.5), f(0.2)],
        [f(-1.0), f(5.0), f(0.8), f(-0.3)],
        [f(0.3), f(-0.7), f(3.0), f(0.6)],
        [f(0.5), f(0.1), f(-0.4), f(6.0)],
    ]);
    let exact = a.solve(Vector(b));

    // Voller Krylov-Raum und Neustart kleiner als die Systemgröße.
    for restart in [4, 2, 1] {
        let mut x = [T::ZERO; 4];
        let mut work = [T::ZERO; 64];
        assert!(gmres_work_len(4, restart).is_some_and(|needed| needed <= work.len()));
        assert!(gmres(&a, &b, &mut x, restart, &mut work, &tol).is_ok());
        assert!(exact.is_ok_and(|e| near(&x, &e.0)));
        let mut image = [T::ZERO; 4];
        assert!(a.apply(&x, &mut image).is_ok() && near(&image, &b));

        // Die Lösung als Startwert: sofort fertig.
        assert!(gmres(&a, &b, &mut x, restart, &mut work, &tol).is_ok_and(|s| s.iterations == 0));
    }

    // Fehlerfälle.
    let mut x = [T::ZERO; 4];
    let mut work = [T::ZERO; 64];
    assert_eq!(gmres(&a, &b, &mut x, 0, &mut work, &tol), Err(SolveError::InvalidInput));
    let mut small = [T::ZERO; 8];
    assert_eq!(gmres(&a, &b, &mut x, 4, &mut small, &tol), Err(SolveError::InvalidInput));
    let mut short = [T::ZERO; 3];
    assert_eq!(gmres(&a, &b, &mut short, 4, &mut work, &tol), Err(SolveError::InvalidInput));
    let none = tol.with_max_iterations(0);
    assert!(matches!(
        gmres(&a, &b, &mut x, 4, &mut work, &none),
        Err(SolveError::NotConverged(s)) if s.iterations == 0
    ));
    assert!(gmres_work_len(usize::MAX, 2).is_none());
}

fn check_preconditioned<T: RealElementary>() {
    let f = T::from_f64;
    let tol = Tolerance::<T>::default().with_rel_x(T::EPSILON.sqrt()).with_max_iterations(50);
    // Symmetrisch positiv definit, aber schlecht skaliert: Die Diagonale
    // überstreicht sechs Größenordnungen.
    let a = Matrix([
        [f(1.0e3), f(1.0), f(0.0), f(0.01)],
        [f(1.0), f(1.0), f(0.01), f(0.0)],
        [f(0.0), f(0.01), f(1.0e-2), f(1.0e-3)],
        [f(0.01), f(0.0), f(1.0e-3), f(1.0e-3)],
    ]);
    let b = [f(1.0), f(-2.0), f(0.5), f(0.25)];
    // Jacobi: Kehrwerte der Diagonale.
    let mut jacobi = Matrix::<T, 4, 4>::ZERO;
    for index in 0..4 {
        jacobi[(index, index)] = a[(index, index)].recip();
    }
    let mut x = [T::ZERO; 4];
    let mut work = [T::ZERO; 16];
    let solved = preconditioned_conjugate_gradient(&a, &jacobi, &b, &mut x, &mut work, &tol);
    assert!(solved.is_ok());
    assert!(a.solve(Vector(b)).is_ok_and(|e| near(&x, &e.0)));
    // Mit der Einheitsmatrix als Vorkonditionierer ist es das gewöhnliche
    // Verfahren und braucht hier mehr Schritte.
    let mut plain = [T::ZERO; 4];
    let unpreconditioned =
        preconditioned_conjugate_gradient(&a, &Matrix::<T, 4, 4>::identity(), &b, &mut plain, &mut work, &tol);
    let steps = |result: &Result<crate::solve::Solution<(), T>, SolveError<(), T>>| match result {
        Ok(state) | Err(SolveError::NotConverged(state)) => state.iterations,
        Err(_) => u32::MAX,
    };
    assert!(steps(&solved) <= steps(&unpreconditioned));

    // Fehlerfälle: zu wenig Arbeitsspeicher, unpassender oder nicht positiv
    // definiter Vorkonditionierer, keine Iterationen erlaubt.
    let mut small = [T::ZERO; 15];
    assert_eq!(
        preconditioned_conjugate_gradient(&a, &jacobi, &b, &mut x, &mut small, &tol),
        Err(SolveError::InvalidInput)
    );
    assert_eq!(
        preconditioned_conjugate_gradient(&a, &Matrix::<T, 3, 3>::identity(), &b, &mut x, &mut work, &tol),
        Err(SolveError::InvalidInput)
    );
    let mut x = [T::ZERO; 4];
    assert!(matches!(
        preconditioned_conjugate_gradient(&a, &(-jacobi), &b, &mut x, &mut work, &tol),
        Err(SolveError::Singular(_))
    ));
    let none = tol.with_max_iterations(0);
    assert!(matches!(
        preconditioned_conjugate_gradient(&a, &jacobi, &b, &mut x, &mut work, &none),
        Err(SolveError::NotConverged(s)) if s.iterations == 0
    ));
    // Null als rechte Seite ist sofort gelöst.
    let mut zero = [T::ZERO; 4];
    assert!(
        preconditioned_conjugate_gradient(&a, &jacobi, &[T::ZERO; 4], &mut zero, &mut work, &tol)
            .is_ok_and(|state| state.iterations == 0)
    );
}

#[test]
fn preconditioned_f32() {
    check_preconditioned::<f32>();
}

#[test]
fn preconditioned_f64() {
    check_preconditioned::<f64>();
}

#[test]
fn gmres_f32() {
    check_gmres::<f32>();
}

#[test]
fn gmres_f64() {
    check_gmres::<f64>();
}

#[test]
fn iterative_f32() {
    check::<f32>();
}

#[test]
fn iterative_f64() {
    check::<f64>();
}
