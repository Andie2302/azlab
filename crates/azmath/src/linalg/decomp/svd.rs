//! Singulärwertzerlegung mit dem einseitigen Jacobi-Verfahren.

use crate::linalg::{LinalgError, Matrix, Vector};
use crate::scalar::RealElementary;

/// Höchstzahl der Durchläufe.
const MAX_SWEEPS: usize = 100;

/// Singulärwertzerlegung `A = U · diag(s) · Vᵀ` einer `R × C`-Matrix mit
/// `R >= C`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Svd<T, const R: usize, const C: usize> {
    /// Linke Singulärvektoren als Spalten (orthonormal, soweit der
    /// Singulärwert nicht null ist).
    pub u: Matrix<T, R, C>,
    /// Singulärwerte, absteigend, nicht negativ.
    pub s: Vector<T, C>,
    /// Rechte Singulärvektoren als Spalten (orthonormal).
    pub v: Matrix<T, C, C>,
}

impl<T: RealElementary, const R: usize, const C: usize> Svd<T, R, C> {
    /// Schwelle, unter der ein Singulärwert als null gilt:
    /// `max(R, C) · EPSILON · größter Singulärwert`.
    fn threshold(&self) -> T {
        let largest = self.s.iter().fold(T::ZERO, |m, &value| m.max(value));
        T::from_f64(R.max(C) as f64) * T::EPSILON * largest
    }

    /// Numerischer Rang: Zahl der Singulärwerte über der Rundungsschwelle.
    pub fn rank(&self) -> usize {
        let threshold = self.threshold();
        self.s.iter().filter(|&&value| value > threshold).count()
    }

    /// Konditionszahl: größter durch kleinsten Singulärwert (unendlich bei
    /// Rangabfall).
    pub fn condition_number(&self) -> T {
        match (self.s.0.first(), self.s.0.last()) {
            (Some(&largest), Some(&smallest)) if smallest > T::ZERO => largest / smallest,
            (Some(_), Some(_)) => T::INFINITY,
            _ => T::ONE,
        }
    }

    /// Pseudoinverse `A⁺ = V · diag(1/s) · Uᵀ`; Singulärwerte unter der
    /// Rundungsschwelle zählen als null.
    pub fn pseudo_inverse(&self) -> Matrix<T, C, R> {
        let threshold = self.threshold();
        let inverse = self.s.map(|value| if value > threshold { value.recip() } else { T::ZERO });
        Matrix::from_fn(|i, j| {
            (0..C).fold(T::ZERO, |sum, k| sum + self.v[(i, k)] * inverse[k] * self.u[(j, k)])
        })
    }

    /// Lösung kleinster Norm von `A·x ≈ b` im Sinn der kleinsten Quadrate;
    /// funktioniert auch bei Rangabfall.
    pub fn solve(&self, b: Vector<T, R>) -> Vector<T, C> {
        self.pseudo_inverse() * b
    }
}

impl<T: RealElementary, const R: usize, const C: usize> Matrix<T, R, C> {
    /// Singulärwertzerlegung mit dem einseitigen Jacobi-Verfahren.
    ///
    /// Braucht `R >= C` (sonst die transponierte Matrix zerlegen). Sehr
    /// genau, aber für kleine Matrizen gedacht.
    // Die Drehungen brauchen den Index für mehrere Felder zugleich.
    #[allow(clippy::needless_range_loop)]
    pub fn svd(self) -> Result<Svd<T, R, C>, LinalgError> {
        if R < C {
            return Err(LinalgError::DimensionMismatch);
        }
        let mut v = Matrix::<T, C, C>::identity().0;
        if !self.norm().is_finite() {
            return Err(LinalgError::NotConverged);
        }
        // Auf Einträge der Größenordnung eins bringen: Die Quadratsummen
        // und ihr Produkt laufen sonst über oder unter und täuschen
        // orthogonale Spalten vor.
        let scale = self.0.iter().flatten().fold(T::ZERO, |largest, &entry| largest.max(entry.abs()));
        if scale == T::ZERO {
            return Ok(Svd { u: self, s: Vector::ZERO, v: Matrix(v) });
        }
        let mut u = self.map(|entry| entry / scale).0;

        let mut converged = false;
        for _ in 0..MAX_SWEEPS {
            let mut rotated = false;
            for p in 0..C {
                for q in p + 1..C {
                    let (mut alpha, mut beta, mut gamma) = (T::ZERO, T::ZERO, T::ZERO);
                    for row in u.iter() {
                        alpha += row[p] * row[p];
                        beta += row[q] * row[q];
                        gamma += row[p] * row[q];
                    }
                    if gamma.abs() <= T::EPSILON * (alpha * beta).sqrt() {
                        continue;
                    }
                    rotated = true;
                    // Drehung, die die Spalten p und q orthogonal macht.
                    let zeta = (beta - alpha) / (gamma + gamma);
                    let t = T::ONE.copysign(zeta) / (zeta.abs() + zeta.hypot(T::ONE));
                    let c = T::ONE.hypot(t).recip();
                    let s = c * t;
                    for row in u.iter_mut() {
                        let (left, right) = (row[p], row[q]);
                        row[p] = c * left - s * right;
                        row[q] = s * left + c * right;
                    }
                    for row in v.iter_mut() {
                        let (left, right) = (row[p], row[q]);
                        row[p] = c * left - s * right;
                        row[q] = s * left + c * right;
                    }
                }
            }
            if !rotated {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(LinalgError::NotConverged);
        }

        // Singulärwerte sind die Längen der Spalten; Spalten normieren.
        let mut s = [T::ZERO; C];
        for col in 0..C {
            let length = Vector::<T, R>(core::array::from_fn(|row| u[row][col])).norm();
            s[col] = length * scale;
            if length > T::ZERO {
                for row in u.iter_mut() {
                    row[col] /= length;
                }
            }
        }
        // Absteigend sortieren.
        for sorted in 0..C {
            let mut largest = sorted;
            for candidate in sorted + 1..C {
                if s[candidate] > s[largest] {
                    largest = candidate;
                }
            }
            if largest != sorted {
                s.swap(sorted, largest);
                for row in u.iter_mut() {
                    row.swap(sorted, largest);
                }
                for row in v.iter_mut() {
                    row.swap(sorted, largest);
                }
            }
        }
        Ok(Svd { u: Matrix(u), s: Vector(s), v: Matrix(v) })
    }
}
