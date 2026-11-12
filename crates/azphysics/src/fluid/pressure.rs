//! Die Druckgleichung: Operator und Vorkonditionierer.

use alloc::vec;
use alloc::vec::Vec;
use azmath::linalg::{LinalgError, LinearOperator};
use azmath::scalar::{Real, RealElementary};

/// Die Zelle hat links eine Nachbarzelle mit Flüssigkeit.
const LEFT: u8 = 1;
/// … rechts.
const RIGHT: u8 = 2;
/// … unten.
const BELOW: u8 = 4;
/// … oben.
const ABOVE: u8 = 8;
/// Die Zelle ist ein Hindernis.
const SOLID: u8 = 16;

/// Aufbau der Druckgleichung für ein Gitter mit Hindernissen.
///
/// Die Matrix ist der negative Laplace-Operator auf den Zellen der
/// Flüssigkeit: Wände und Hindernisse lassen nichts durch (kein Eintrag),
/// jenseits offener Ränder ist der Druck null, Hindernis-Zellen bilden auf
/// sich selbst ab. Sie ist symmetrisch und positiv semidefinit; mit
/// mindestens einem offenen Rand ist sie positiv definit.
#[derive(Clone, Debug)]
pub(super) struct Stencil<T> {
    nx: usize,
    /// Je Zelle: welche Nachbarn Flüssigkeit enthalten.
    mask: Vec<u8>,
    /// Diagonale der Matrix: Zahl der durchlässigen Seiten.
    diagonal: Vec<T>,
    /// Kehrwerte der Diagonale des unvollständigen Cholesky-Faktors.
    factors: Vec<T>,
}

impl<T: RealElementary> Stencil<T> {
    /// Speicher für `nx × ny` Zellen; vor Gebrauch [`rebuild`](Self::rebuild)
    /// aufrufen.
    pub(super) fn new(nx: usize, ny: usize) -> Self {
        let cells = nx * ny;
        Self { nx, mask: vec![SOLID; cells], diagonal: vec![T::ONE; cells], factors: vec![T::ZERO; cells] }
    }

    /// Baut Matrix und Vorkonditionierer aus den Hindernissen und den
    /// offenen Rändern (links, rechts, unten, oben).
    ///
    /// Der Vorkonditionierer ist die modifizierte unvollständige
    /// Cholesky-Zerlegung ohne Auffüllung, MIC(0), nach R. Bridson, „Fluid
    /// Simulation for Computer Graphics“, Abschnitt 5.5.
    pub(super) fn rebuild(&mut self, solid: &[bool], open: [bool; 4]) {
        let nx = self.nx;
        let cells = self.mask.len();
        let ny = cells / nx.max(1);
        for index in 0..cells {
            if solid[index] {
                self.mask[index] = SOLID;
                self.diagonal[index] = T::ONE;
                continue;
            }
            let (i, j) = (index % nx, index / nx);
            let mut mask = 0;
            let mut sides = 0.0;
            for (inside, offset, bit, open) in [
                (i > 0, index.wrapping_sub(1), LEFT, open[0]),
                (i + 1 < nx, index + 1, RIGHT, open[1]),
                (j > 0, index.wrapping_sub(nx), BELOW, open[2]),
                (j + 1 < ny, index + nx, ABOVE, open[3]),
            ] {
                if inside && !solid[offset] {
                    mask |= bit;
                    sides += 1.0;
                } else if !inside && open {
                    sides += 1.0;
                }
            }
            self.mask[index] = mask;
            self.diagonal[index] = T::from_f64(sides);
        }

        // Anteil der weggelassenen Einträge, der auf die Diagonale kommt,
        // und Schutz gegen zu kleine Diagonalwerte.
        let (tuning, safety) = (T::from_f64(0.97), T::from_f64(0.25));
        for index in 0..cells {
            let mask = self.mask[index];
            if mask & SOLID != 0 {
                self.factors[index] = T::ZERO;
                continue;
            }
            let mut pivot = self.diagonal[index];
            if mask & LEFT != 0 {
                let factor = self.factors[index - 1];
                pivot -= factor * factor;
                // Der linke Nachbar würde mit seinem oberen einen Eintrag
                // erzeugen, den es in der Matrix nicht gibt.
                if self.mask[index - 1] & ABOVE != 0 {
                    pivot -= tuning * factor * factor;
                }
            }
            if mask & BELOW != 0 {
                let factor = self.factors[index - nx];
                pivot -= factor * factor;
                if self.mask[index - nx] & RIGHT != 0 {
                    pivot -= tuning * factor * factor;
                }
            }
            if pivot < safety * self.diagonal[index] {
                pivot = self.diagonal[index];
            }
            self.factors[index] = if pivot > T::ZERO { pivot.sqrt().recip() } else { T::ZERO };
        }
    }
}

/// Die Matrix der Druckgleichung.
pub(super) struct PressureOperator<'a, T>(pub(super) &'a Stencil<T>);

impl<T: Real> LinearOperator<T> for PressureOperator<'_, T> {
    fn rows(&self) -> usize {
        self.0.mask.len()
    }

    fn cols(&self) -> usize {
        self.0.mask.len()
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        let Stencil { nx, mask, diagonal, .. } = self.0;
        if x.len() != mask.len() || y.len() != mask.len() {
            return Err(LinalgError::DimensionMismatch);
        }
        for (index, out) in y.iter_mut().enumerate() {
            let mask = mask[index];
            if mask & SOLID != 0 {
                *out = x[index];
                continue;
            }
            let mut sum = diagonal[index] * x[index];
            if mask & LEFT != 0 {
                sum -= x[index - 1];
            }
            if mask & RIGHT != 0 {
                sum -= x[index + 1];
            }
            if mask & BELOW != 0 {
                sum -= x[index - nx];
            }
            if mask & ABOVE != 0 {
                sum -= x[index + nx];
            }
            *out = sum;
        }
        Ok(())
    }
}

/// Vorkonditionierer zur Druckgleichung: wendet die Inverse der
/// unvollständigen Cholesky-Zerlegung an (vorwärts und rückwärts
/// einsetzen).
pub(super) struct Preconditioner<'a, T>(pub(super) &'a Stencil<T>);

impl<T: Real> LinearOperator<T> for Preconditioner<'_, T> {
    fn rows(&self) -> usize {
        self.0.mask.len()
    }

    fn cols(&self) -> usize {
        self.0.mask.len()
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        let Stencil { nx, mask, factors, .. } = self.0;
        let size = mask.len();
        if x.len() != size || y.len() != size {
            return Err(LinalgError::DimensionMismatch);
        }
        // Vorwärts: L·q = x.
        for index in 0..size {
            let mask = mask[index];
            let mut sum = x[index];
            if mask & SOLID != 0 {
                y[index] = sum;
                continue;
            }
            if mask & LEFT != 0 {
                sum += factors[index - 1] * y[index - 1];
            }
            if mask & BELOW != 0 {
                sum += factors[index - nx] * y[index - nx];
            }
            y[index] = sum * factors[index];
        }
        // Rückwärts: Lᵀ·z = q.
        for index in (0..size).rev() {
            let mask = mask[index];
            if mask & SOLID != 0 {
                continue;
            }
            let mut sum = y[index];
            if mask & RIGHT != 0 {
                sum += factors[index] * y[index + 1];
            }
            if mask & ABOVE != 0 {
                sum += factors[index] * y[index + nx];
            }
            y[index] = sum * factors[index];
        }
        Ok(())
    }
}
