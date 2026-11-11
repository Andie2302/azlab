//! Besitzender Tensor (Feature `alloc`).

use super::ops;
use super::view::{TensorView, TensorViewMut, element_count};
use crate::linalg::LinalgError;
use crate::scalar::Algebra;
use alloc::vec;
use alloc::vec::Vec;

/// Tensor der Stufe `RANK` auf dem Heap; die Form steht erst zur Laufzeit
/// fest.
#[derive(Clone, Debug, PartialEq)]
pub struct DTensor<T, const RANK: usize> {
    data: Vec<T>,
    shape: [usize; RANK],
}

impl<T: Algebra, const RANK: usize> DTensor<T, RANK> {
    /// Nulltensor der Form `shape`. Fehler, wenn die Zahl der Elemente nicht
    /// darstellbar ist.
    pub fn zeros(shape: [usize; RANK]) -> Result<Self, LinalgError> {
        let len = element_count(&shape).ok_or(LinalgError::DimensionMismatch)?;
        Ok(Self { data: vec![T::ZERO; len], shape })
    }

    /// Tensor aus hintereinander abgelegten Elementen (letzter Index läuft am
    /// schnellsten).
    ///
    /// Fehler, wenn `data` nicht zur Form passt.
    pub fn from_vec(data: Vec<T>, shape: [usize; RANK]) -> Result<Self, LinalgError> {
        if element_count(&shape) == Some(data.len()) {
            Ok(Self { data, shape })
        } else {
            Err(LinalgError::DimensionMismatch)
        }
    }

    /// Die Form: Ausdehnung je Index.
    #[inline]
    pub fn shape(&self) -> [usize; RANK] {
        self.shape
    }

    /// Lesende Sicht.
    #[inline]
    pub fn view(&self) -> TensorView<'_, T, RANK> {
        TensorView::from_parts(&self.data, self.shape)
    }

    /// Schreibende Sicht.
    #[inline]
    pub fn view_mut(&mut self) -> TensorViewMut<'_, T, RANK> {
        TensorViewMut::from_parts(&mut self.data, self.shape)
    }

    /// Gibt die Elemente als `Vec` heraus.
    #[inline]
    pub fn into_vec(self) -> Vec<T> {
        self.data
    }

    /// Tensor aus den Ausdehnungen, die `dims` liefert. Fehler, wenn es
    /// nicht genau `RANK` sind.
    fn with_dims(dims: impl Iterator<Item = usize>) -> Result<Self, LinalgError> {
        let mut shape = [0_usize; RANK];
        let mut count = 0_usize;
        for dim in dims {
            let Some(slot) = shape.get_mut(count) else {
                return Err(LinalgError::DimensionMismatch);
            };
            *slot = dim;
            count += 1;
        }
        if count != RANK {
            return Err(LinalgError::DimensionMismatch);
        }
        Self::zeros(shape)
    }

    /// Äußeres Produkt, siehe [`outer`](super::outer). Die Stufe des
    /// Ergebnisses muss die Summe der Stufen der Eingaben sein.
    pub fn outer<const RA: usize, const RB: usize>(
        a: TensorView<'_, T, RA>,
        b: TensorView<'_, T, RB>,
    ) -> Result<Self, LinalgError> {
        let mut out = Self::with_dims(a.shape().into_iter().chain(b.shape()))?;
        ops::outer(a, b, &mut out.data)?;
        Ok(out)
    }

    /// Verallgemeinertes Produkt, siehe [`dot`](super::dot). Die Stufe des
    /// Ergebnisses muss die Summe der Stufen der Eingaben minus zwei sein.
    pub fn dot<const RA: usize, const RB: usize>(
        a: TensorView<'_, T, RA>,
        b: TensorView<'_, T, RB>,
    ) -> Result<Self, LinalgError> {
        let leading = a.shape().into_iter().take(RA.saturating_sub(1));
        let trailing = b.shape().into_iter().skip(1);
        let mut out = Self::with_dims(leading.chain(trailing))?;
        ops::dot(a, b, &mut out.data)?;
        Ok(out)
    }

    /// Kontraktion, siehe [`contract`](super::contract). Die Stufe des
    /// Ergebnisses muss um zwei kleiner sein als die der Eingabe.
    pub fn contract<const RA: usize>(
        a: TensorView<'_, T, RA>,
        first: usize,
        second: usize,
    ) -> Result<Self, LinalgError> {
        let remaining = a
            .shape()
            .into_iter()
            .enumerate()
            .filter(|&(axis, _)| axis != first && axis != second)
            .map(|(_, dim)| dim);
        let mut out = Self::with_dims(remaining)?;
        ops::contract(a, first, second, &mut out.data)?;
        Ok(out)
    }
}
