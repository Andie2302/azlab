//! Finite-Differenzen-Operatoren auf regelmäßigen Gittern.

use super::grid::UniformGrid;
use crate::linalg::{LinalgError, LinearOperator};
use crate::scalar::Real;

/// Randbedingung einer Raumrichtung: Was gilt für den Nachbarn jenseits des
/// Gitterrands?
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Boundary {
    /// Wert null außerhalb. Das Gitter enthält dann nur die inneren Punkte;
    /// der Rand liegt einen Gitterabstand jenseits des ersten und letzten.
    Dirichlet,
    /// Kein Fluss über den Rand: Der Nachbar außerhalb hat den Wert des
    /// Randpunkts.
    Neumann,
    /// Periodisch: Hinter dem letzten Punkt folgt wieder der erste.
    Periodic,
}

/// Differenzenformel für die erste Ableitung.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scheme {
    /// Zentral, 2. Ordnung: `(rechts − links) / (2h)`.
    Central,
    /// Vorwärts, 1. Ordnung: `(rechts − mitte) / h`. Aufwind-Formel für
    /// Transport nach links (negative Geschwindigkeit).
    Forward,
    /// Rückwärts, 1. Ordnung: `(mitte − links) / h`. Aufwind-Formel für
    /// Transport nach rechts (positive Geschwindigkeit).
    Backward,
}

/// Form eines Gitters mit den Schrittweiten im Puffer je Richtung.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Layout<const D: usize> {
    shape: [usize; D],
    strides: [usize; D],
    len: usize,
}

impl<const D: usize> Layout<D> {
    fn new(shape: [usize; D]) -> Self {
        let mut strides = [1_usize; D];
        let mut len = 1_usize;
        for (stride, &extent) in strides.iter_mut().zip(&shape).rev() {
            *stride = len;
            len = len.saturating_mul(extent);
        }
        Self { shape, strides, len }
    }

    /// Werte links und rechts von Position `flat` in Richtung `axis`.
    fn neighbours<T: Real>(&self, u: &[T], flat: usize, axis: usize, boundary: Boundary) -> (T, T) {
        let (stride, extent) = (self.strides[axis], self.shape[axis]);
        let coordinate = (flat / stride) % extent;
        let centre = u[flat];
        let outside = |wrapped: usize| match boundary {
            Boundary::Dirichlet => T::ZERO,
            Boundary::Neumann => centre,
            Boundary::Periodic => u[wrapped],
        };
        let left = if coordinate > 0 {
            u[flat - stride]
        } else {
            outside(flat + stride * (extent - 1))
        };
        let right = if coordinate + 1 < extent {
            u[flat + stride]
        } else {
            outside(flat - stride * (extent - 1))
        };
        (left, right)
    }
}

/// Laplace-Operator `Δ = Σ ∂²/∂x_d²` mit dem Standard-Differenzenstern
/// 2. Ordnung.
///
/// Mit Dirichlet-Rändern ist `−Δ` symmetrisch positiv definit und damit für
/// [`conjugate_gradient`](crate::solve::conjugate_gradient) geeignet; `−Δ`
/// erhält man über
/// [`AffineOperator`](crate::linalg::AffineOperator)`::new(&laplace, 0, −1)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Laplacian<T, const D: usize> {
    layout: Layout<D>,
    inverse_square: [T; D],
    boundary: [Boundary; D],
}

impl<T: Real, const D: usize> Laplacian<T, D> {
    /// Laplace-Operator auf `grid` mit je einer Randbedingung pro Richtung.
    pub fn new(grid: &UniformGrid<T, D>, boundary: [Boundary; D]) -> Self {
        Self {
            layout: Layout::new(grid.shape()),
            inverse_square: grid.spacing().map(|h| (h * h).recip()),
            boundary,
        }
    }
}

impl<T: Real, const D: usize> LinearOperator<T> for Laplacian<T, D> {
    #[inline]
    fn rows(&self) -> usize {
        self.layout.len
    }

    #[inline]
    fn cols(&self) -> usize {
        self.layout.len
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if x.len() != self.layout.len || y.len() != self.layout.len {
            return Err(LinalgError::DimensionMismatch);
        }
        let two = T::from_f64(2.0);
        for (flat, out) in y.iter_mut().enumerate() {
            let centre = x[flat];
            let mut sum = T::ZERO;
            for axis in 0..D {
                let (left, right) = self.layout.neighbours(x, flat, axis, self.boundary[axis]);
                sum += (left - two * centre + right) * self.inverse_square[axis];
            }
            *out = sum;
        }
        Ok(())
    }
}

/// Erste Ableitung `∂/∂x_axis` als Differenzenoperator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Derivative<T, const D: usize> {
    layout: Layout<D>,
    axis: usize,
    inverse_spacing: T,
    scheme: Scheme,
    boundary: Boundary,
}

impl<T: Real, const D: usize> Derivative<T, D> {
    /// Ableitung in Richtung `axis` auf `grid`.
    ///
    /// `None`, wenn `axis >= D`.
    pub fn new(
        grid: &UniformGrid<T, D>,
        axis: usize,
        scheme: Scheme,
        boundary: Boundary,
    ) -> Option<Self> {
        let spacing = *grid.spacing().get(axis)?;
        Some(Self {
            layout: Layout::new(grid.shape()),
            axis,
            inverse_spacing: spacing.recip(),
            scheme,
            boundary,
        })
    }
}

impl<T: Real, const D: usize> LinearOperator<T> for Derivative<T, D> {
    #[inline]
    fn rows(&self) -> usize {
        self.layout.len
    }

    #[inline]
    fn cols(&self) -> usize {
        self.layout.len
    }

    fn apply(&self, x: &[T], y: &mut [T]) -> Result<(), LinalgError> {
        if x.len() != self.layout.len || y.len() != self.layout.len {
            return Err(LinalgError::DimensionMismatch);
        }
        for (flat, out) in y.iter_mut().enumerate() {
            let (left, right) = self.layout.neighbours(x, flat, self.axis, self.boundary);
            let difference = match self.scheme {
                Scheme::Central => (right - left) * T::from_f64(0.5),
                Scheme::Forward => right - x[flat],
                Scheme::Backward => x[flat] - left,
            };
            *out = difference * self.inverse_spacing;
        }
        Ok(())
    }
}
