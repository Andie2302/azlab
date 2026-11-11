//! Eigenwerte: symmetrische Matrizen mit Eigenvektoren (Jacobi-Verfahren),
//! allgemeine Matrizen nur Eigenwerte (QR-Verfahren).
//!
//! Beide Verfahren sind für kleine Matrizen fester Größe gedacht.

use crate::linalg::{LinalgError, Matrix, Vector};
use crate::number::Complex;
use crate::poly::solve_quadratic_complex;
use crate::scalar::{Algebra, Conjugate, RealElementary};

/// Höchstzahl der Durchläufe des Jacobi-Verfahrens.
const MAX_SWEEPS: usize = 100;

/// Eigenwerte und Eigenvektoren einer symmetrischen Matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SymmetricEigen<T, const N: usize> {
    /// Eigenwerte, aufsteigend.
    pub values: Vector<T, N>,
    /// Eigenvektoren als Spalten, orthonormal; Spalte `j` gehört zu
    /// `values[j]`.
    pub vectors: Matrix<T, N, N>,
}

impl<T: RealElementary, const N: usize> Matrix<T, N, N> {
    /// Eigenwerte und Eigenvektoren einer symmetrischen Matrix mit dem
    /// zyklischen Jacobi-Verfahren.
    ///
    /// Ob die Matrix symmetrisch ist, wird nicht geprüft. Die Eigenvektoren
    /// sind bis auf Rundung orthonormal; es gilt
    /// `A = V · diag(values) · Vᵀ`.
    // Die Drehungen brauchen den Index für Zeilen und Spalten zugleich.
    #[allow(clippy::needless_range_loop)]
    pub fn symmetric_eigen(self) -> Result<SymmetricEigen<T, N>, LinalgError> {
        let mut a = self.0;
        let mut v = Matrix::<T, N, N>::identity().0;
        let scale = self.norm();
        if !scale.is_finite() {
            return Err(LinalgError::NotConverged);
        }
        let threshold = T::EPSILON * scale;

        let mut converged = false;
        for _ in 0..MAX_SWEEPS {
            let mut off = T::ZERO;
            for p in 0..N {
                for q in p + 1..N {
                    off = off.max(a[p][q].abs());
                }
            }
            if off <= threshold {
                converged = true;
                break;
            }
            for p in 0..N {
                for q in p + 1..N {
                    if a[p][q].abs() <= threshold * T::EPSILON {
                        continue;
                    }
                    // Drehwinkel, der a[p][q] auslöscht.
                    let theta = (a[q][q] - a[p][p]) / (a[p][q] + a[p][q]);
                    let t = T::ONE.copysign(theta) / (theta.abs() + theta.hypot(T::ONE));
                    let c = T::ONE.hypot(t).recip();
                    let s = t * c;
                    for k in 0..N {
                        let (left, right) = (a[k][p], a[k][q]);
                        a[k][p] = c * left - s * right;
                        a[k][q] = s * left + c * right;
                    }
                    for k in 0..N {
                        let (upper, lower) = (a[p][k], a[q][k]);
                        a[p][k] = c * upper - s * lower;
                        a[q][k] = s * upper + c * lower;
                    }
                    for k in 0..N {
                        let (left, right) = (v[k][p], v[k][q]);
                        v[k][p] = c * left - s * right;
                        v[k][q] = s * left + c * right;
                    }
                }
            }
        }
        if !converged {
            return Err(LinalgError::NotConverged);
        }

        // Aufsteigend sortieren, Spalten der Eigenvektoren mitnehmen.
        let mut values: [T; N] = core::array::from_fn(|i| a[i][i]);
        for sorted in 0..N {
            let mut smallest = sorted;
            for candidate in sorted + 1..N {
                if values[candidate] < values[smallest] {
                    smallest = candidate;
                }
            }
            if smallest != sorted {
                values.swap(sorted, smallest);
                for row in v.iter_mut() {
                    row.swap(sorted, smallest);
                }
            }
        }
        Ok(SymmetricEigen { values: Vector(values), vectors: Matrix(v) })
    }

    /// Alle Eigenwerte einer beliebigen reellen Matrix.
    ///
    /// Sie können komplex sein und treten dann in konjugierten Paaren auf.
    /// Die Reihenfolge ist nicht festgelegt. Für symmetrische Matrizen ist
    /// [`symmetric_eigen`](Self::symmetric_eigen) genauer und liefert auch
    /// die Eigenvektoren.
    pub fn eigenvalues(self) -> Result<[Complex<T>; N], LinalgError> {
        self.map(Complex::from_real).eigenvalues()
    }
}

impl<T: RealElementary, const N: usize> Matrix<Complex<T>, N, N> {
    /// Alle Eigenwerte einer komplexen Matrix mit dem QR-Verfahren
    /// (Hessenberg-Form, Givens-Drehungen, Wilkinson-Verschiebung).
    ///
    /// Die Reihenfolge ist nicht festgelegt.
    // Die Drehungen brauchen den Index für Zeilen und Spalten zugleich.
    #[allow(clippy::needless_range_loop)]
    pub fn eigenvalues(self) -> Result<[Complex<T>; N], LinalgError> {
        let zero = Complex::<T>::ZERO;
        let one = Complex::<T>::ONE;
        let mut a = self.0;
        let scale = self.norm();
        if !scale.is_finite() {
            return Err(LinalgError::NotConverged);
        }

        // Drehung, die (x, y) auf (r, 0) abbildet: liefert (c, s) mit
        // conj(c)·x + conj(s)·y = r und −s·x + c·y = 0.
        let rotation = |x: Complex<T>, y: Complex<T>| {
            let radius = x.norm().hypot(y.norm());
            if radius == T::ZERO { (one, zero) } else { (x / radius, y / radius) }
        };

        // Auf Hessenberg-Form bringen (Ähnlichkeitstransformation).
        for col in 0..N.saturating_sub(2) {
            for row in col + 2..N {
                let (c, s) = rotation(a[col + 1][col], a[row][col]);
                for k in 0..N {
                    let (upper, lower) = (a[col + 1][k], a[row][k]);
                    a[col + 1][k] = c.conj() * upper + s.conj() * lower;
                    a[row][k] = c * lower - s * upper;
                }
                for k in 0..N {
                    let (left, right) = (a[k][col + 1], a[k][row]);
                    a[k][col + 1] = left * c + right * s;
                    a[k][row] = right * c.conj() - left * s.conj();
                }
            }
        }

        // QR-Schritte auf dem aktiven Block [0, active); unten rechts fertige
        // Eigenwerte werden abgespalten.
        let mut active = N;
        let mut stalled = 0_usize;
        let mut budget = 100_usize.saturating_mul(N.max(1));
        let mut rotations = [(one, zero); N];
        while active > 1 {
            let (last, before) = (active - 1, active - 2);
            let size = a[last][last].abs_sum() + a[before][before].abs_sum();
            let reference = if size > T::ZERO { size } else { scale };
            if a[last][before].abs_sum() <= T::EPSILON * reference {
                a[last][before] = zero;
                active -= 1;
                stalled = 0;
                continue;
            }
            if budget == 0 {
                return Err(LinalgError::NotConverged);
            }
            budget -= 1;
            stalled += 1;

            // Wilkinson-Verschiebung: der Eigenwert des letzten 2×2-Blocks,
            // der näher am letzten Diagonalelement liegt.
            let trace = a[before][before] + a[last][last];
            let determinant = a[before][before] * a[last][last] - a[before][last] * a[last][before];
            let [first, second] = solve_quadratic_complex([determinant, -trace, one]);
            let corner = a[last][last];
            let mut shift = if (first - corner).abs_sum() <= (second - corner).abs_sum() {
                first
            } else {
                second
            };
            if stalled % 10 == 0 {
                // Ausnahmeverschiebung gegen Stillstand.
                shift += Complex::from_real(a[last][before].abs_sum());
            }

            for k in 0..active {
                a[k][k] -= shift;
            }
            // R = Gₘ…G₁·A
            for k in 0..last {
                let (c, s) = rotation(a[k][k], a[k + 1][k]);
                rotations[k] = (c, s);
                for col in k..active {
                    let (upper, lower) = (a[k][col], a[k + 1][col]);
                    a[k][col] = c.conj() * upper + s.conj() * lower;
                    a[k + 1][col] = c * lower - s * upper;
                }
            }
            // A ← R·G₁ᴴ…Gₘᴴ
            for k in 0..last {
                let (c, s) = rotations[k];
                for row in 0..(k + 2).min(active) {
                    let (left, right) = (a[row][k], a[row][k + 1]);
                    a[row][k] = left * c + right * s;
                    a[row][k + 1] = right * c.conj() - left * s.conj();
                }
            }
            for k in 0..active {
                a[k][k] += shift;
            }
        }
        Ok(core::array::from_fn(|i| a[i][i]))
    }
}
