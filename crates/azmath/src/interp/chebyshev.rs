//! Näherung von Funktionen durch Tschebyschow-Reihen.

use crate::poly::clenshaw;
use crate::scalar::RealElementary;

/// Die `N` Tschebyschow-Stützstellen auf [a, b], aufsteigend.
///
/// Interpolation an diesen Stellen vermeidet das Aufschwingen, das bei
/// gleichabständigen Stützstellen auftritt.
pub fn chebyshev_nodes<T: RealElementary, const N: usize>(a: T, b: T) -> [T; N] {
    let half = T::from_f64(0.5);
    let (mid, radius) = ((a + b) * half, (b - a) * half);
    let count = T::from_f64(N as f64);
    core::array::from_fn(|index| {
        // Von rechts nach links nummeriert, damit das Ergebnis aufsteigt.
        let k = T::from_f64((N - 1 - index) as f64);
        mid + radius * (T::PI * (k + half) / count).cos()
    })
}

/// Tschebyschow-Reihe mit `N` Gliedern auf dem Intervall [a, b].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChebyshevSeries<T, const N: usize> {
    coefficients: [T; N],
    a: T,
    b: T,
}

impl<T: RealElementary, const N: usize> ChebyshevSeries<T, N> {
    /// Nähert `f` auf [a, b] durch die Reihe an, die `f` an den
    /// Tschebyschow-Stützstellen interpoliert.
    ///
    /// Für glatte Funktionen fällt der Fehler mit wachsendem `N` sehr
    /// schnell.
    pub fn fit(mut f: impl FnMut(T) -> T, a: T, b: T) -> Self {
        let half = T::from_f64(0.5);
        let (mid, radius) = ((a + b) * half, (b - a) * half);
        let count = T::from_f64(N as f64);
        let angle = |k: usize| T::PI * (T::from_f64(k as f64) + half) / count;
        let samples: [T; N] = core::array::from_fn(|k| f(mid + radius * angle(k).cos()));
        let coefficients = core::array::from_fn(|j| {
            let order = T::from_f64(j as f64);
            let sum = samples
                .iter()
                .enumerate()
                .fold(T::ZERO, |acc, (k, &sample)| acc + sample * (order * angle(k)).cos());
            let scale = if j == 0 { T::ONE } else { T::from_f64(2.0) };
            sum * scale / count
        });
        Self { coefficients, a, b }
    }

    /// Die Koeffizienten `c[j]` der Reihe `Σ c[j]·T_j`.
    #[inline]
    pub fn coefficients(&self) -> &[T; N] {
        &self.coefficients
    }

    /// Wert der Reihe an der Stelle `x` (gedacht für `x` in [a, b]).
    pub fn eval(&self, x: T) -> T {
        let two = T::from_f64(2.0);
        let scaled = (two * x - self.a - self.b) / (self.b - self.a);
        clenshaw(&self.coefficients, scaled)
    }
}
