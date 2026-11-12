//! Transport entlang der Strömung (Semi-Lagrange, wahlweise mit
//! MacCormack-Korrektur).

use alloc::vec;
use alloc::vec::Vec;
use azmath::scalar::RealElementary;

/// Lage der Werte eines Felds im Gitter: `width × height` Werte, der Wert
/// `(i, j)` sitzt bei `(i + offset_x, j + offset_y)` in Zellbreiten.
#[derive(Clone, Copy, Debug)]
pub(super) struct Layout {
    pub(super) width: usize,
    pub(super) height: usize,
    pub(super) offset_x: f64,
    pub(super) offset_y: f64,
}

/// Zerlegt eine Koordinate in den Index des linken Stützpunkts und den
/// Anteil zum rechten; außerhalb wird auf den Rand begrenzt.
#[inline]
fn split<T: RealElementary>(coordinate: T, offset: f64, count: usize) -> (usize, usize, T) {
    let last = count.saturating_sub(1);
    let position = coordinate.to_f64() - offset;
    // Auch NaN landet so bei null.
    let clamped = if position > 0.0 { position.min(last as f64) } else { 0.0 };
    // Abschneiden statt `floor`: Der Wert ist nicht negativ.
    let low = (clamped as usize).min(last.saturating_sub(1));
    (low, (low + 1).min(last), T::from_f64(clamped - low as f64))
}

impl Layout {
    /// Horizontale Geschwindigkeit: auf den senkrechten Zellkanten.
    pub(super) fn u(nx: usize, ny: usize) -> Self {
        Self { width: nx + 1, height: ny, offset_x: 0.0, offset_y: 0.5 }
    }

    /// Vertikale Geschwindigkeit: auf den waagrechten Zellkanten.
    pub(super) fn v(nx: usize, ny: usize) -> Self {
        Self { width: nx, height: ny + 1, offset_x: 0.5, offset_y: 0.0 }
    }

    /// Größen in der Zellmitte.
    pub(super) fn centre(nx: usize, ny: usize) -> Self {
        Self { width: nx, height: ny, offset_x: 0.5, offset_y: 0.5 }
    }

    /// Bilinear interpolierter Wert an der Stelle `(x, y)` (in Zellbreiten)
    /// samt kleinstem und größtem der vier beteiligten Werte.
    pub(super) fn sample<T: RealElementary>(&self, data: &[T], x: T, y: T) -> (T, T, T) {
        let (i0, i1, fx) = split(x, self.offset_x, self.width);
        let (j0, j1, fy) = split(y, self.offset_y, self.height);
        let a = data[j0 * self.width + i0];
        let b = data[j0 * self.width + i1];
        let c = data[j1 * self.width + i0];
        let d = data[j1 * self.width + i1];
        let bottom = a + (b - a) * fx;
        let top = c + (d - c) * fx;
        (bottom + (top - bottom) * fy, a.min(b).min(c.min(d)), a.max(b).max(c.max(d)))
    }
}

/// Das Geschwindigkeitsfeld, entlang dessen transportiert wird.
pub(super) struct Flow<'a, T> {
    pub(super) u: &'a [T],
    pub(super) v: &'a [T],
    pub(super) nx: usize,
    pub(super) ny: usize,
    /// Zeitschritt geteilt durch die Zellbreite.
    pub(super) step: T,
}

impl<T: RealElementary> Flow<'_, T> {
    /// Geschwindigkeit an der Stelle `(x, y)`.
    fn at(&self, x: T, y: T) -> (T, T) {
        (Layout::u(self.nx, self.ny).sample(self.u, x, y).0, Layout::v(self.nx, self.ny).sample(self.v, x, y).0)
    }

    /// Verfolgt die Strömung von `(x, y)` aus um einen Zeitschritt;
    /// `direction` ist −1 für rückwärts, +1 für vorwärts. Mit `midpoint`
    /// nach der Mittelpunktsregel, sonst mit einem Euler-Schritt.
    fn trace(&self, x: T, y: T, direction: T, midpoint: bool) -> (T, T) {
        let step = self.step * direction;
        let (mut u, mut v) = self.at(x, y);
        if midpoint {
            let half = step * T::from_f64(0.5);
            (u, v) = self.at(x + half * u, y + half * v);
        }
        (x + step * u, y + step * v)
    }

    /// Transportiert `source` nach `target`.
    ///
    /// Semi-Lagrange: Jeder Gitterpunkt holt sich den Wert von dort, wo die
    /// Strömung ihn herträgt. Das ist für jede Schrittweite stabil, aber
    /// verschmiert. Mit `maccormack` wird der Fehler aus einem Schritt hin
    /// und zurück geschätzt und abgezogen (Selle u. a. 2008); die Werte
    /// bleiben dabei zwischen denen der Nachbarn am Herkunftsort. Die
    /// Korrektur gleicht auch den Fehler der Bahnverfolgung aus, deshalb
    /// genügt dort ein Euler-Schritt.
    pub(super) fn advect(
        &self,
        layout: Layout,
        source: &[T],
        target: &mut [T],
        scratch: &mut Scratch<T>,
        maccormack: bool,
    ) {
        let position = |index: usize| {
            (
                T::from_f64((index % layout.width) as f64 + layout.offset_x),
                T::from_f64((index / layout.width) as f64 + layout.offset_y),
            )
        };
        if !maccormack {
            for (index, out) in target.iter_mut().enumerate().take(source.len()) {
                let (x, y) = position(index);
                let (bx, by) = self.trace(x, y, -T::ONE, true);
                *out = layout.sample(source, bx, by).0;
            }
            return;
        }
        for index in 0..source.len() {
            let (x, y) = position(index);
            let (bx, by) = self.trace(x, y, -T::ONE, false);
            (scratch.first[index], scratch.low[index], scratch.high[index]) = layout.sample(source, bx, by);
        }
        let half = T::from_f64(0.5);
        for (index, out) in target.iter_mut().enumerate().take(source.len()) {
            let (x, y) = position(index);
            let (fx, fy) = self.trace(x, y, T::ONE, false);
            let returned = layout.sample(&scratch.first, fx, fy).0;
            let corrected = scratch.first[index] + (source[index] - returned) * half;
            *out = corrected.max(scratch.low[index]).min(scratch.high[index]);
        }
    }
}

/// Arbeitsspeicher für den Transport.
#[derive(Clone, Debug)]
pub(super) struct Scratch<T> {
    /// Ergebnis des einfachen Transports.
    first: Vec<T>,
    /// Kleinster und größter Wert am Herkunftsort.
    low: Vec<T>,
    high: Vec<T>,
}

impl<T: RealElementary> Scratch<T> {
    /// Speicher für Felder bis zur Länge `len`.
    pub(super) fn new(len: usize) -> Self {
        Self { first: vec![T::ZERO; len], low: vec![T::ZERO; len], high: vec![T::ZERO; len] }
    }
}
