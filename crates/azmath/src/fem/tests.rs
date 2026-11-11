//! Tests der Finite-Elemente-Bausteine.

use super::{
    ElementOperator, FemError, PlaneModel, assemble_load, line2_load, line2_mass,
    line2_stiffness, tetra4_load, tetra4_mass, tetra4_stiffness, triangle3_elasticity,
    triangle3_load, triangle3_mass, triangle3_stiffness,
};
use crate::linalg::{LinalgError, LinearOperator, Matrix, Vector};
use crate::scalar::Real;
use crate::testutil::close;

fn close_m<T: Real, const N: usize>(a: Matrix<T, N, N>, b: Matrix<T, N, N>, ulps: f64) -> bool {
    (0..N).all(|i| (0..N).all(|j| close(a[(i, j)], b[(i, j)], ulps)))
}

/// Summe aller Einträge.
fn total<T: Real, const N: usize>(m: Matrix<T, N, N>) -> T {
    m.0.iter().flatten().fold(T::ZERO, |sum, &value| sum + value)
}

/// Zeilensummen null (Konstanten liegen im Kern) und symmetrisch.
fn is_stiffness<T: Real, const N: usize>(m: Matrix<T, N, N>, ulps: f64) -> bool {
    (m * Vector::splat(T::ONE)).iter().all(|&row| close(row, T::ZERO, ulps))
        && close_m(m, m.transpose(), ulps)
}

fn elements<T: Real>() {
    let f = T::from_f64;

    // Strecke der Länge 2.
    let line = [[f(1.0)], [f(3.0)]];
    assert_eq!(line2_stiffness(&line), Ok(Matrix([[f(0.5), f(-0.5)], [f(-0.5), f(0.5)]])));
    assert!(line2_mass(&line).is_ok_and(|m| close(total(m), f(2.0), 8.0) && m[(0, 0)] == m[(0, 1)] * f(2.0)));
    assert_eq!(line2_load(&line), Ok(Vector([f(1.0), f(1.0)])));
    assert_eq!(line2_stiffness(&[[f(1.0)], [f(1.0)]]), Err(FemError::DegenerateElement));

    // Rechtwinkliges Einheitsdreieck: bekannte Steifigkeitsmatrix.
    let triangle = [[f(0.0), f(0.0)], [f(1.0), f(0.0)], [f(0.0), f(1.0)]];
    assert_eq!(
        triangle3_stiffness(&triangle),
        Ok(Matrix([
            [f(1.0), f(-0.5), f(-0.5)],
            [f(-0.5), f(0.5), f(0.0)],
            [f(-0.5), f(0.0), f(0.5)],
        ]))
    );
    assert!(triangle3_mass(&triangle).is_ok_and(|m| close(total(m), f(0.5), 8.0)));
    assert!(triangle3_load(&triangle).is_ok_and(|l| close(l.sum(), f(0.5), 8.0)));
    // Beliebiges Dreieck, beide Umlaufrichtungen.
    let skew = [[f(0.2), f(-0.1)], [f(1.7), f(0.4)], [f(0.6), f(1.3)]];
    let flipped = [skew[0], skew[2], skew[1]];
    assert!(triangle3_stiffness(&skew).is_ok_and(|k| is_stiffness(k, 64.0)));
    assert!(match (triangle3_mass(&skew), triangle3_mass(&flipped)) {
        (Ok(a), Ok(b)) => close(total(a), total(b), 16.0),
        _ => false,
    });
    let collinear = [[f(0.0), f(0.0)], [f(1.0), f(1.0)], [f(2.0), f(2.0)]];
    assert_eq!(triangle3_stiffness(&collinear), Err(FemError::DegenerateElement));
    assert_eq!(triangle3_mass(&collinear), Err(FemError::DegenerateElement));
    assert_eq!(triangle3_load(&collinear), Err(FemError::DegenerateElement));

    // Einheitstetraeder: Volumen 1/6.
    let tetra = [
        [f(0.0), f(0.0), f(0.0)],
        [f(1.0), f(0.0), f(0.0)],
        [f(0.0), f(1.0), f(0.0)],
        [f(0.0), f(0.0), f(1.0)],
    ];
    let sixth = T::ONE / f(6.0);
    assert!(tetra4_stiffness(&tetra).is_ok_and(|k| {
        is_stiffness(k, 64.0) && close(k[(0, 0)], f(0.5), 16.0) && close(k[(1, 1)], sixth, 16.0)
    }));
    assert!(tetra4_mass(&tetra).is_ok_and(|m| close(total(m), sixth, 16.0)));
    assert!(tetra4_load(&tetra).is_ok_and(|l| close(l.sum(), sixth, 16.0)));
    let flat = [tetra[0], tetra[1], tetra[2], [f(1.0), f(1.0), f(0.0)]];
    assert_eq!(tetra4_stiffness(&flat), Err(FemError::DegenerateElement));

    // Ebene Elastizität: symmetrisch, und starre Bewegungen (zwei
    // Verschiebungen, eine kleine Drehung) erzeugen keine Kräfte.
    for model in [PlaneModel::Stress, PlaneModel::Strain] {
        let stiffness = triangle3_elasticity(&skew, f(210.0), f(0.3), f(0.5), model);
        assert!(stiffness.is_ok_and(|k| {
            let translate_x = Vector([T::ONE, T::ZERO, T::ONE, T::ZERO, T::ONE, T::ZERO]);
            let translate_y = Vector([T::ZERO, T::ONE, T::ZERO, T::ONE, T::ZERO, T::ONE]);
            let rotate = Vector([
                -skew[0][1], skew[0][0], -skew[1][1], skew[1][0], -skew[2][1], skew[2][0],
            ]);
            let stretch = Vector([skew[0][0], T::ZERO, skew[1][0], T::ZERO, skew[2][0], T::ZERO]);
            close_m(k, k.transpose(), 4096.0)
                && [translate_x, translate_y, rotate]
                    .into_iter()
                    .all(|mode| (k * mode).iter().all(|&force| force.abs() < f(1.0e-2)))
                // Eine Dehnung kostet dagegen Energie.
                && stretch.dot(k * stretch) > T::ONE
        }));
    }
    assert_eq!(
        triangle3_elasticity(&collinear, f(1.0), f(0.3), f(1.0), PlaneModel::Stress),
        Err(FemError::DegenerateElement)
    );
}

fn operators<T: Real>() {
    let f = T::from_f64;
    // Stab aus vier Elementen auf [0, 1].
    let nodes = [[f(0.0)], [f(0.25)], [f(0.5)], [f(0.75)], [f(1.0)]];
    let cells = [[0, 1], [1, 2], [2, 3], [3, 4]];

    // Ohne Randbedingung: Tridiagonalmatrix (1/h)·[−1, 2, −1].
    let free = ElementOperator::new(&nodes, &cells, None, line2_stiffness);
    assert!(free.is_ok_and(|operator| {
        let x = [f(1.0), f(3.0), f(-2.0), f(0.5), f(4.0)];
        let mut y = [T::ZERO; 5];
        let four = f(4.0);
        operator.rows() == 5
            && operator.cols() == 5
            && operator.apply(&x, &mut y).is_ok()
            && y[0] == four * (x[0] - x[1])
            && y[2] == four * (f(2.0) * x[2] - x[1] - x[3])
            && y[4] == four * (x[4] - x[3])
            && operator.apply(&x[..4], &mut y) == Err(LinalgError::DimensionMismatch)
    }));

    // Mit festgehaltenen Enden: dort Einheitsmatrix, keine Kopplung.
    let fixed = [true, false, false, false, true];
    let constrained = ElementOperator::new(&nodes, &cells, Some(&fixed), line2_stiffness);
    assert!(constrained.as_ref().is_ok_and(|operator| {
        let x = [f(7.0), f(3.0), f(-2.0), f(0.5), f(9.0)];
        let mut y = [T::ZERO; 5];
        let four = f(4.0);
        operator.apply(&x, &mut y).is_ok()
            && y[0] == x[0]
            && y[4] == x[4]
            && y[1] == four * (f(2.0) * x[1] - x[2])
            && y[3] == four * (f(2.0) * x[3] - x[2])
    }));
    // Vorgaben eintragen: u(0) = 1, u(1) = 3 ohne Quelle → Gerade.
    assert!(constrained.is_ok_and(|operator| {
        let values = [f(1.0), T::ZERO, T::ZERO, T::ZERO, f(3.0)];
        let mut rhs = [T::ZERO; 5];
        let mut work = [T::ZERO; 5];
        let line = [f(1.0), f(1.5), f(2.0), f(2.5), f(3.0)];
        let mut image = [T::ZERO; 5];
        operator.lift(&values, &mut rhs, &mut work).is_ok()
            && operator.apply(&line, &mut image).is_ok()
            && image.iter().zip(rhs).all(|(&got, want)| close(got, want, 16.0))
            && operator.lift(&values[..4], &mut rhs, &mut work) == Err(FemError::LengthMismatch)
    }));

    // Rechte Seite: Jeder innere Knoten bekommt h, jeder Randknoten h/2.
    let mut load = [T::ZERO; 5];
    assert!(assemble_load(&nodes, &cells, line2_load, &mut load).is_ok());
    assert_eq!(load, [f(0.125), f(0.25), f(0.25), f(0.25), f(0.125)]);
    assert_eq!(
        assemble_load(&nodes, &cells, line2_load, &mut load[..4]),
        Err(FemError::LengthMismatch)
    );

    // Fehlerhafte Netze.
    let broken = [[0, 9]];
    assert!(matches!(
        ElementOperator::new(&nodes, &broken, None, line2_stiffness),
        Err(FemError::InvalidMesh)
    ));
    let folded = [[1, 1]];
    assert!(matches!(
        ElementOperator::new(&nodes, &folded, None, line2_stiffness),
        Err(FemError::DegenerateElement)
    ));
    assert!(matches!(
        ElementOperator::new(&nodes, &cells, Some(&fixed[..3]), line2_stiffness),
        Err(FemError::LengthMismatch)
    ));
    assert_eq!(assemble_load(&nodes, &broken, line2_load, &mut load), Err(FemError::InvalidMesh));

    // Zwei Freiheitsgrade je Knoten: Elastizität auf einem Dreieck.
    let corners = [[f(0.0), f(0.0)], [f(1.0), f(0.0)], [f(0.0), f(1.0)]];
    let single = [[0, 1, 2]];
    let elastic = |c: &[[T; 2]; 3]| triangle3_elasticity(c, f(1.0), f(0.25), f(1.0), PlaneModel::Stress);
    assert!(ElementOperator::new(&corners, &single, None, elastic).is_ok_and(|operator| {
        let x = [f(0.1), f(-0.2), f(0.3), f(0.0), f(-0.1), f(0.2)];
        let mut y = [T::ZERO; 6];
        operator.rows() == 6
            && operator.apply(&x, &mut y).is_ok()
            && elastic(&corners).is_ok_and(|k| y == (k * Vector(x)).0)
    }));
}

#[test]
fn fem_f32() {
    elements::<f32>();
    operators::<f32>();
}

#[test]
fn fem_f64() {
    elements::<f64>();
    operators::<f64>();
}

#[cfg(any(feature = "std", feature = "libm"))]
mod problems {
    use super::*;
    use crate::fem::{quad4_mass, quad4_stiffness};
    use crate::scalar::RealElementary;
    use crate::solve::{Tolerance, conjugate_gradient};

    const CELLS: usize = 8;
    const NODES: usize = (CELLS + 1) * (CELLS + 1);

    fn quads<T: RealElementary>() {
        let f = T::from_f64;
        // Einheitsquadrat: bekannte Werte 2/3 auf der Diagonale, −1/3 gegenüber.
        let square = [[f(0.0), f(0.0)], [f(1.0), f(0.0)], [f(1.0), f(1.0)], [f(0.0), f(1.0)]];
        assert!(quad4_stiffness(&square).is_ok_and(|k| {
            is_stiffness(k, 64.0)
                && close(k[(0, 0)], f(2.0) / f(3.0), 16.0)
                && close(k[(0, 2)], -T::ONE / f(3.0), 16.0)
                && close(k[(0, 1)], -T::ONE / f(6.0), 16.0)
        }));
        assert!(quad4_mass(&square).is_ok_and(|m| {
            close(total(m), T::ONE, 16.0) && close(m[(0, 0)], T::ONE / f(9.0), 16.0)
        }));
        // Allgemeines Viereck: Masse = Fläche (hier ein Trapez mit Fläche 3).
        let trapezoid = [[f(0.0), f(0.0)], [f(2.0), f(0.0)], [f(2.0), f(2.0)], [f(0.0), f(1.0)]];
        assert!(quad4_mass(&trapezoid).is_ok_and(|m| close(total(m), f(3.0), 16.0)));
        assert!(quad4_stiffness(&trapezoid).is_ok_and(|k| is_stiffness(k, 64.0)));
        // Im Uhrzeigersinn nummeriert: verdreht.
        let clockwise = [square[0], square[3], square[2], square[1]];
        assert_eq!(quad4_stiffness(&clockwise), Err(FemError::DegenerateElement));
    }

    /// −Δu = 2π² sin(πx)·sin(πy) auf dem Einheitsquadrat mit u = 0 am Rand;
    /// Lösung sin(πx)·sin(πy). Dreiecksnetz aus 8 × 8 Zellen.
    fn poisson<T: RealElementary>() {
        let f = T::from_f64;
        let pi = T::PI;
        let step = T::ONE / f(CELLS as f64);
        let nodes: [[T; 2]; NODES] = core::array::from_fn(|index| {
            [step * f((index % (CELLS + 1)) as f64), step * f((index / (CELLS + 1)) as f64)]
        });
        // Jede Zelle in zwei Dreiecke teilen.
        let triangles: [[usize; 3]; 2 * CELLS * CELLS] = core::array::from_fn(|index| {
            let (cell, upper) = (index / 2, index % 2 == 1);
            let corner = (cell / CELLS) * (CELLS + 1) + cell % CELLS;
            if upper {
                [corner + 1, corner + CELLS + 2, corner + CELLS + 1]
            } else {
                [corner, corner + 1, corner + CELLS + 1]
            }
        });
        let fixed: [bool; NODES] = core::array::from_fn(|index| {
            let (col, row) = (index % (CELLS + 1), index / (CELLS + 1));
            col == 0 || row == 0 || col == CELLS || row == CELLS
        });
        let source = |x: T, y: T| f(2.0) * pi * pi * (pi * x).sin() * (pi * y).sin();

        // Rechte Seite: Quelle am Schwerpunkt jedes Elements.
        let mut rhs = [T::ZERO; NODES];
        let loaded = assemble_load(
            &nodes,
            &triangles,
            |corners: &[[T; 2]; 3]| {
                let third = T::ONE / f(3.0);
                let x = (corners[0][0] + corners[1][0] + corners[2][0]) * third;
                let y = (corners[0][1] + corners[1][1] + corners[2][1]) * third;
                triangle3_load(corners).map(|load| load * source(x, y))
            },
            &mut rhs,
        );
        assert!(loaded.is_ok());

        let operator = ElementOperator::new(&nodes, &triangles, Some(&fixed), triangle3_stiffness);
        let mut u = [T::ZERO; NODES];
        let solved = operator.is_ok_and(|operator| {
            let mut work = [T::ZERO; 3 * NODES];
            let tol = Tolerance::<T>::default()
                .with_rel_x(T::EPSILON.sqrt())
                .with_max_iterations(500);
            let boundary_values = [T::ZERO; NODES];
            let mut scratch = [T::ZERO; NODES];
            operator.lift(&boundary_values, &mut rhs, &mut scratch).is_ok()
                && conjugate_gradient(&operator, &rhs, &mut u, &mut work, &tol).is_ok()
        });
        assert!(solved);
        let worst = nodes.iter().zip(u).fold(T::ZERO, |worst, (&[x, y], value)| {
            worst.max((value - (pi * x).sin() * (pi * y).sin()).abs())
        });
        assert!(worst < f(0.06));
        // Randwerte bleiben null.
        assert!(fixed.iter().zip(u).all(|(&is_fixed, value)| !is_fixed || value == T::ZERO));
    }

    #[test]
    fn problems_f32() {
        quads::<f32>();
        poisson::<f32>();
    }

    #[test]
    fn problems_f64() {
        quads::<f64>();
        poisson::<f64>();
    }
}

#[cfg(feature = "alloc")]
mod assembled {
    use super::*;
    use crate::fem::assemble_csr;

    fn check<T: Real>() {
        let f = T::from_f64;
        let nodes = [[f(0.0), f(0.0)], [f(1.0), f(0.0)], [f(1.0), f(1.0)], [f(0.0), f(1.0)], [f(0.5), f(0.5)]];
        let triangles = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]];
        let x = [f(1.0), f(-2.0), f(0.5), f(3.0), f(-1.5)];

        let assembled = assemble_csr(&nodes, &triangles, triangle3_stiffness);
        let operator = ElementOperator::new(&nodes, &triangles, None, triangle3_stiffness);
        let ok = match (assembled, operator) {
            (Ok(csr), Ok(operator)) => {
                let (mut from_csr, mut from_elements) = ([T::ZERO; 5], [T::ZERO; 5]);
                csr.view().is_ok_and(|view| {
                    view.apply(&x, &mut from_csr).is_ok()
                        // Mittelknoten ist mit allen verbunden, Ecken nicht mit der gegenüber.
                        && view.get(0, 2) == Some(T::ZERO)
                        && view.get(4, 4).is_some_and(|d| d > T::ZERO)
                }) && operator.apply(&x, &mut from_elements).is_ok()
                    && from_csr.iter().zip(from_elements).all(|(&a, b)| close(a, b, 64.0))
                    // 5 Diagonal-, 8 Rand- und 8 Mittelkopplungen.
                    && csr.nnz() == 21
                    && csr.clone() == csr
            }
            _ => false,
        };
        assert!(ok);
        assert!(matches!(
            assemble_csr(&nodes, &[[0, 1, 7]], triangle3_stiffness),
            Err(FemError::InvalidMesh)
        ));
        assert!(matches!(
            assemble_csr(&nodes, &[[0, 0, 1]], triangle3_stiffness),
            Err(FemError::DegenerateElement)
        ));
        // Knoten ohne Element ergeben leere Zeilen.
        let lonely = [[f(0.0)], [f(1.0)], [f(5.0)]];
        assert!(assemble_csr(&lonely, &[[0, 1]], line2_mass)
            .is_ok_and(|csr| csr.nnz() == 4 && csr.view().is_ok_and(|v| v.get(2, 2) == Some(T::ZERO))));
    }

    #[test]
    fn assembled_f32() {
        check::<f32>();
    }

    #[test]
    fn assembled_f64() {
        check::<f64>();
    }
}
