//! Tensor-Sichten auf Puffer des Aufrufers.

use crate::linalg::LinalgError;
use core::ops::{Index, IndexMut};

/// Anzahl der Elemente zu einer Form; `None` bei Überlauf.
pub(super) fn element_count<const RANK: usize>(shape: &[usize; RANK]) -> Option<usize> {
    shape.iter().try_fold(1_usize, |count, &dim| count.checked_mul(dim))
}

/// Position im Puffer; `None`, wenn ein Index außerhalb liegt.
fn offset<const RANK: usize>(shape: &[usize; RANK], index: [usize; RANK]) -> Option<usize> {
    let mut position = 0_usize;
    for (&dim, component) in shape.iter().zip(index) {
        if component >= dim {
            return None;
        }
        position = position.checked_mul(dim)?.checked_add(component)?;
    }
    Some(position)
}

/// Lesende Sicht auf einen Tensor der Stufe `RANK`.
#[derive(Clone, Copy, Debug)]
pub struct TensorView<'a, T, const RANK: usize> {
    data: &'a [T],
    shape: [usize; RANK],
}

/// Schreibende Sicht auf einen Tensor der Stufe `RANK`.
#[derive(Debug)]
pub struct TensorViewMut<'a, T, const RANK: usize> {
    data: &'a mut [T],
    shape: [usize; RANK],
}

impl<'a, T, const RANK: usize> TensorView<'a, T, RANK> {
    /// Sicht auf `data` als Tensor der Form `shape`.
    ///
    /// Fehler, wenn `data` nicht genau so viele Elemente hat wie das Produkt
    /// der Form.
    pub fn new(data: &'a [T], shape: [usize; RANK]) -> Result<Self, LinalgError> {
        if element_count(&shape) == Some(data.len()) {
            Ok(Self { data, shape })
        } else {
            Err(LinalgError::DimensionMismatch)
        }
    }

    /// Sicht ohne Längenprüfung; der Aufrufer garantiert die passende Länge.
    #[cfg(feature = "alloc")]
    #[inline]
    pub(super) fn from_parts(data: &'a [T], shape: [usize; RANK]) -> Self {
        Self { data, shape }
    }

    /// Die Form: Ausdehnung je Index.
    #[inline]
    pub fn shape(&self) -> [usize; RANK] {
        self.shape
    }

    /// Anzahl der Elemente.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// `true`, wenn der Tensor keine Elemente hat.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Die Elemente hintereinander; der letzte Index läuft am schnellsten.
    #[inline]
    pub fn data(&self) -> &'a [T] {
        self.data
    }

    /// Geprüfter Zugriff: `None`, wenn ein Index außerhalb liegt.
    #[inline]
    pub fn get(&self, index: [usize; RANK]) -> Option<&T> {
        self.data.get(offset(&self.shape, index)?)
    }

    /// Dieselben Daten mit anderer Form, auch anderer Stufe.
    ///
    /// Fehler, wenn die neue Form eine andere Zahl von Elementen hat.
    pub fn reshape<const NEW: usize>(
        self,
        shape: [usize; NEW],
    ) -> Result<TensorView<'a, T, NEW>, LinalgError> {
        TensorView::new(self.data, shape)
    }
}

impl<'a, T, const RANK: usize> TensorViewMut<'a, T, RANK> {
    /// Sicht auf `data` als Tensor der Form `shape`.
    ///
    /// Fehler, wenn `data` nicht genau so viele Elemente hat wie das Produkt
    /// der Form.
    pub fn new(data: &'a mut [T], shape: [usize; RANK]) -> Result<Self, LinalgError> {
        if element_count(&shape) == Some(data.len()) {
            Ok(Self { data, shape })
        } else {
            Err(LinalgError::DimensionMismatch)
        }
    }

    /// Sicht ohne Längenprüfung; der Aufrufer garantiert die passende Länge.
    #[cfg(feature = "alloc")]
    #[inline]
    pub(super) fn from_parts(data: &'a mut [T], shape: [usize; RANK]) -> Self {
        Self { data, shape }
    }

    /// Die Form: Ausdehnung je Index.
    #[inline]
    pub fn shape(&self) -> [usize; RANK] {
        self.shape
    }

    /// Anzahl der Elemente.
    #[inline]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// `true`, wenn der Tensor keine Elemente hat.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Geprüfter Zugriff: `None`, wenn ein Index außerhalb liegt.
    #[inline]
    pub fn get(&self, index: [usize; RANK]) -> Option<&T> {
        self.data.get(offset(&self.shape, index)?)
    }

    /// Geprüfter schreibender Zugriff.
    #[inline]
    pub fn get_mut(&mut self, index: [usize; RANK]) -> Option<&mut T> {
        self.data.get_mut(offset(&self.shape, index)?)
    }

    /// Lesende Sicht auf dieselben Daten.
    #[inline]
    pub fn as_view(&self) -> TensorView<'_, T, RANK> {
        TensorView { data: self.data, shape: self.shape }
    }
}

impl<T: Copy, const RANK: usize> TensorViewMut<'_, T, RANK> {
    /// Ersetzt jedes Element `x` durch `f(x)`.
    pub fn map_in_place(&mut self, mut f: impl FnMut(T) -> T) {
        for slot in self.data.iter_mut() {
            *slot = f(*slot);
        }
    }

    /// Ersetzt jedes Element `x` durch `f(x, y)` mit dem entsprechenden
    /// Element `y` von `other`.
    ///
    /// Fehler, wenn die Formen verschieden sind.
    pub fn zip_in_place(
        &mut self,
        other: TensorView<'_, T, RANK>,
        mut f: impl FnMut(T, T) -> T,
    ) -> Result<(), LinalgError> {
        if self.shape != other.shape {
            return Err(LinalgError::DimensionMismatch);
        }
        for (slot, &rhs) in self.data.iter_mut().zip(other.data) {
            *slot = f(*slot, rhs);
        }
        Ok(())
    }
}

/// Zugriff `t[[i, j, k]]` wie bei Arrays: bricht bei falschem Index ab
/// (Panic). Für geprüften Zugriff gibt es [`TensorView::get`].
impl<T, const RANK: usize> Index<[usize; RANK]> for TensorView<'_, T, RANK> {
    type Output = T;

    #[inline]
    fn index(&self, index: [usize; RANK]) -> &T {
        &self.data[offset(&self.shape, index).unwrap_or(usize::MAX)]
    }
}

/// Zugriff `t[[i, j, k]]` wie bei Arrays: bricht bei falschem Index ab
/// (Panic). Für geprüften Zugriff gibt es [`TensorViewMut::get`].
impl<T, const RANK: usize> Index<[usize; RANK]> for TensorViewMut<'_, T, RANK> {
    type Output = T;

    #[inline]
    fn index(&self, index: [usize; RANK]) -> &T {
        &self.data[offset(&self.shape, index).unwrap_or(usize::MAX)]
    }
}

/// Schreibender Zugriff wie bei Arrays: bricht bei falschem Index ab (Panic).
/// Für geprüften Zugriff gibt es [`TensorViewMut::get_mut`].
impl<T, const RANK: usize> IndexMut<[usize; RANK]> for TensorViewMut<'_, T, RANK> {
    #[inline]
    fn index_mut(&mut self, index: [usize; RANK]) -> &mut T {
        &mut self.data[offset(&self.shape, index).unwrap_or(usize::MAX)]
    }
}
