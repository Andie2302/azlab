//! Regelmäßige Gitter.

use crate::scalar::Real;

/// Regelmäßiges Gitter in `D` Raumrichtungen: `shape[d]` Punkte je Richtung
/// im Abstand `spacing[d]`, beginnend bei `origin`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformGrid<T, const D: usize> {
    shape: [usize; D],
    origin: [T; D],
    spacing: [T; D],
}

impl<T: Real, const D: usize> UniformGrid<T, D> {
    /// Erzeugt das Gitter.
    ///
    /// `None`, wenn die Gesamtzahl der Punkte nicht darstellbar ist oder ein
    /// Abstand nicht positiv und endlich ist.
    pub fn new(shape: [usize; D], origin: [T; D], spacing: [T; D]) -> Option<Self> {
        shape.iter().try_fold(1_usize, |count, &extent| count.checked_mul(extent))?;
        if spacing.iter().all(|&h| h > T::ZERO && h.is_finite()) {
            Some(Self { shape, origin, spacing })
        } else {
            None
        }
    }

    /// Anzahl der Punkte je Richtung.
    #[inline]
    pub fn shape(&self) -> [usize; D] {
        self.shape
    }

    /// Gitterabstand je Richtung.
    #[inline]
    pub fn spacing(&self) -> [T; D] {
        self.spacing
    }

    /// Gesamtzahl der Gitterpunkte.
    #[inline]
    pub fn len(&self) -> usize {
        self.shape.iter().product()
    }

    /// `true`, wenn das Gitter keine Punkte hat.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Position des Gitterpunkts `index` im Puffer; `None` außerhalb.
    pub fn flat(&self, index: [usize; D]) -> Option<usize> {
        let mut position = 0_usize;
        for (&extent, component) in self.shape.iter().zip(index) {
            if component >= extent {
                return None;
            }
            position = position * extent + component;
        }
        Some(position)
    }

    /// Gitterindex zur Position `flat` im Puffer; `None` außerhalb.
    pub fn index_of(&self, flat: usize) -> Option<[usize; D]> {
        if flat >= self.len() {
            return None;
        }
        let mut rest = flat;
        let mut index = [0_usize; D];
        for (component, &extent) in index.iter_mut().zip(&self.shape).rev() {
            *component = rest % extent;
            rest /= extent;
        }
        Some(index)
    }

    /// Koordinaten des Gitterpunkts `index`.
    pub fn point(&self, index: [usize; D]) -> [T; D] {
        core::array::from_fn(|d| self.origin[d] + self.spacing[d] * T::from_f64(index[d] as f64))
    }

    /// Füllt `field` mit den Werten von `f` an allen Gitterpunkten.
    ///
    /// Gibt `false` zurück (und ändert nichts), wenn `field` nicht genau
    /// [`len`](Self::len) Elemente hat.
    pub fn sample(&self, field: &mut [T], mut f: impl FnMut([T; D]) -> T) -> bool {
        if field.len() != self.len() {
            return false;
        }
        for (flat, slot) in field.iter_mut().enumerate() {
            if let Some(index) = self.index_of(flat) {
                *slot = f(self.point(index));
            }
        }
        true
    }
}
