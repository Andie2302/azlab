//! Gauß-Legendre-Quadratur.

use crate::scalar::RealElementary;

/// Gauß-Legendre-Regel mit `N` Stützstellen.
///
/// Integriert Polynome bis zum Grad `2N − 1` exakt. Stützstellen und
/// Gewichte werden beim Erzeugen berechnet (Nullstellen des
/// Legendre-Polynoms per Newton-Verfahren), nicht aus Tabellen gelesen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GaussLegendre<T, const N: usize> {
    nodes: [T; N],
    weights: [T; N],
}

impl<T: RealElementary, const N: usize> GaussLegendre<T, N> {
    /// Berechnet Stützstellen und Gewichte auf dem Intervall [−1, 1].
    pub fn new() -> Self {
        let mut nodes = [T::ZERO; N];
        let mut weights = [T::ZERO; N];
        let count = T::from_f64(N as f64);
        // Legendre-Polynom P_N und seine Ableitung über die Rekursion.
        let legendre = |x: T| {
            let (mut current, mut previous) = (T::ONE, T::ZERO);
            for degree in 1..=N {
                let k = T::from_f64(degree as f64);
                let next = ((k + k - T::ONE) * x * current - (k - T::ONE) * previous) / k;
                previous = current;
                current = next;
            }
            (current, count * (x * current - previous) / (x * x - T::ONE))
        };
        for i in 0..N.div_ceil(2) {
            // Startwert: Näherung der i-ten Nullstelle.
            let mut x = (T::PI * (T::from_f64(i as f64) + T::from_f64(0.75))
                / (count + T::from_f64(0.5)))
            .cos();
            for _ in 0..100 {
                let (value, slope) = legendre(x);
                let step = value / slope;
                x -= step;
                if step.abs() <= T::EPSILON * T::from_f64(4.0) {
                    break;
                }
            }
            // Ableitung an der endgültigen Nullstelle für das Gewicht.
            let (_, slope) = legendre(x);
            let weight = T::from_f64(2.0) / ((T::ONE - x * x) * slope * slope);
            nodes[i] = -x;
            nodes[N - 1 - i] = x;
            weights[i] = weight;
            weights[N - 1 - i] = weight;
        }
        Self { nodes, weights }
    }

    /// Stützstellen auf [−1, 1], aufsteigend.
    #[inline]
    pub fn nodes(&self) -> &[T; N] {
        &self.nodes
    }

    /// Gewichte zu den Stützstellen.
    #[inline]
    pub fn weights(&self) -> &[T; N] {
        &self.weights
    }

    /// Integral von `f` über [a, b].
    pub fn integrate(&self, mut f: impl FnMut(T) -> T, a: T, b: T) -> T {
        let half = (b - a) * T::from_f64(0.5);
        let mid = (a + b) * T::from_f64(0.5);
        let sum = self
            .nodes
            .iter()
            .zip(&self.weights)
            .fold(T::ZERO, |acc, (&x, &w)| acc + w * f(mid + half * x));
        sum * half
    }

    /// Integral von `f` über [a, b], aufgeteilt in `panels` gleich breite
    /// Teilintervalle (mindestens eines), auf jedem die Regel angewandt.
    pub fn integrate_panels(&self, mut f: impl FnMut(T) -> T, a: T, b: T, panels: u32) -> T {
        let panels = panels.max(1);
        let width = (b - a) / T::from_f64(f64::from(panels));
        let mut sum = T::ZERO;
        for index in 0..panels {
            let left = a + width * T::from_f64(f64::from(index));
            sum += self.integrate(&mut f, left, left + width);
        }
        sum
    }
}

impl<T: RealElementary, const N: usize> Default for GaussLegendre<T, N> {
    fn default() -> Self {
        Self::new()
    }
}
