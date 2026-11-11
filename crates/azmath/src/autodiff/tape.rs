//! Aufzeichnung der Rechenschritte für den Rückwärtsmodus.

use crate::scalar::Algebra;
use core::cell::RefCell;
use core::fmt;

/// Nummer für „nicht auf dem Band“: Konstanten und fehlende Vorgänger.
pub(super) const DETACHED: usize = usize::MAX;

/// Fehler des Rückwärtsmodus.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum TapeError {
    /// Das Band war zu klein für die Rechnung.
    CapacityExceeded,
    /// Eine Variable gehört nicht zu diesem Band (z. B. nach `clear`).
    UnknownVariable,
}

impl fmt::Display for TapeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::CapacityExceeded => "Band ist zu klein",
            Self::UnknownVariable => "Variable gehört nicht zu diesem Band",
        })
    }
}

impl core::error::Error for TapeError {}

/// Ein Rechenschritt: bis zu zwei Vorgänger und die partiellen Ableitungen
/// des Ergebnisses nach ihnen.
#[derive(Clone, Copy, Debug)]
pub(super) struct Node<T> {
    parents: [usize; 2],
    partials: [T; 2],
}

/// Rückwärtsdurchlauf über die Schritte `nodes`: verteilt die Ableitung des
/// Ausgangs auf alle Vorgänger.
fn sweep<T: Algebra>(nodes: &[Node<T>], adjoints: &mut [T], output: usize) {
    adjoints.fill(T::ZERO);
    if let Some(slot) = adjoints.get_mut(output) {
        *slot = T::ONE;
    }
    for index in (0..nodes.len().min(output.saturating_add(1))).rev() {
        let adjoint = adjoints[index];
        if adjoint == T::ZERO {
            continue;
        }
        let node = nodes[index];
        for (parent, partial) in node.parents.into_iter().zip(node.partials) {
            if let Some(slot) = adjoints.get_mut(parent) {
                *slot += partial * adjoint;
            }
        }
    }
}

/// Speicher für die Rechenschritte des Rückwärtsmodus.
///
/// Implementiert von [`Tape`] (feste Größe, ohne Heap) und – mit dem Feature
/// `alloc` – von `GrowTape`. Anwender benutzen die Methoden der Bänder und
/// von [`Var`](super::Var), nicht diesen Trait.
pub trait Record<T> {
    /// Zeichnet einen Schritt auf und gibt seine Nummer zurück.
    #[doc(hidden)]
    fn push(&self, parents: [usize; 2], partials: [T; 2]) -> usize;

    /// Ableitungen des Schritts `output` nach den Schritten `wanted`.
    #[doc(hidden)]
    fn backward(&self, output: usize, wanted: &[usize], out: &mut [T]) -> Result<(), TapeError>;
}

struct Fixed<T, const CAP: usize> {
    nodes: [Node<T>; CAP],
    adjoints: [T; CAP],
    len: usize,
    overflowed: bool,
}

/// Band fester Größe für bis zu `CAP` Rechenschritte, ohne Heap.
///
/// Jede Variable und jede Rechenoperation belegt einen Platz. Reicht der
/// Platz nicht, meldet die Auswertung der Ableitung
/// [`TapeError::CapacityExceeded`]; die Rechnung selbst läuft mit richtigen
/// Werten weiter.
pub struct Tape<T, const CAP: usize> {
    inner: RefCell<Fixed<T, CAP>>,
}

impl<T: Algebra, const CAP: usize> Tape<T, CAP> {
    /// Leeres Band.
    pub fn new() -> Self {
        let empty = Node { parents: [DETACHED; 2], partials: [T::ZERO; 2] };
        Self {
            inner: RefCell::new(Fixed {
                nodes: [empty; CAP],
                adjoints: [T::ZERO; CAP],
                len: 0,
                overflowed: false,
            }),
        }
    }

    /// Anzahl der belegten Plätze.
    pub fn len(&self) -> usize {
        self.inner.borrow().len
    }

    /// `true`, wenn noch nichts aufgezeichnet ist.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Löscht das Band. Zuvor erzeugte Variablen sind danach ungültig.
    pub fn clear(&self) {
        let mut inner = self.inner.borrow_mut();
        inner.len = 0;
        inner.overflowed = false;
    }
}

impl<T: Algebra, const CAP: usize> Default for Tape<T, CAP> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T, const CAP: usize> fmt::Debug for Tape<T, CAP> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Tape").field("capacity", &CAP).finish_non_exhaustive()
    }
}

impl<T: Algebra, const CAP: usize> Record<T> for Tape<T, CAP> {
    fn push(&self, parents: [usize; 2], partials: [T; 2]) -> usize {
        let mut inner = self.inner.borrow_mut();
        let index = inner.len;
        match inner.nodes.get_mut(index) {
            Some(slot) => {
                *slot = Node { parents, partials };
                inner.len += 1;
                index
            }
            None => {
                inner.overflowed = true;
                DETACHED
            }
        }
    }

    fn backward(&self, output: usize, wanted: &[usize], out: &mut [T]) -> Result<(), TapeError> {
        let mut inner = self.inner.borrow_mut();
        if inner.overflowed {
            return Err(TapeError::CapacityExceeded);
        }
        let len = inner.len;
        if (output != DETACHED && output >= len) || wanted.iter().any(|&w| w != DETACHED && w >= len) {
            return Err(TapeError::UnknownVariable);
        }
        let Fixed { nodes, adjoints, .. } = &mut *inner;
        sweep(&nodes[..len], &mut adjoints[..len], output);
        for (slot, &index) in out.iter_mut().zip(wanted) {
            *slot = adjoints.get(index).copied().unwrap_or(T::ZERO);
        }
        Ok(())
    }
}

#[cfg(feature = "alloc")]
mod growable {
    use super::{DETACHED, Node, Record, TapeError, sweep};
    use crate::scalar::Algebra;
    use alloc::vec;
    use alloc::vec::Vec;
    use core::cell::RefCell;
    use core::fmt;

    /// Band, das mit der Rechnung wächst (Feature `alloc`).
    pub struct GrowTape<T> {
        nodes: RefCell<Vec<Node<T>>>,
    }

    impl<T: Algebra> GrowTape<T> {
        /// Leeres Band.
        pub fn new() -> Self {
            Self { nodes: RefCell::new(Vec::new()) }
        }

        /// Anzahl der aufgezeichneten Schritte.
        pub fn len(&self) -> usize {
            self.nodes.borrow().len()
        }

        /// `true`, wenn noch nichts aufgezeichnet ist.
        pub fn is_empty(&self) -> bool {
            self.len() == 0
        }

        /// Löscht das Band. Zuvor erzeugte Variablen sind danach ungültig.
        pub fn clear(&self) {
            self.nodes.borrow_mut().clear();
        }
    }

    impl<T: Algebra> Default for GrowTape<T> {
        fn default() -> Self {
            Self::new()
        }
    }

    impl<T> fmt::Debug for GrowTape<T> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.debug_struct("GrowTape").finish_non_exhaustive()
        }
    }

    impl<T: Algebra> Record<T> for GrowTape<T> {
        fn push(&self, parents: [usize; 2], partials: [T; 2]) -> usize {
            let mut nodes = self.nodes.borrow_mut();
            nodes.push(Node { parents, partials });
            nodes.len() - 1
        }

        fn backward(&self, output: usize, wanted: &[usize], out: &mut [T]) -> Result<(), TapeError> {
            let nodes = self.nodes.borrow();
            let len = nodes.len();
            if (output != DETACHED && output >= len)
                || wanted.iter().any(|&w| w != DETACHED && w >= len)
            {
                return Err(TapeError::UnknownVariable);
            }
            let mut adjoints = vec![T::ZERO; len];
            sweep(&nodes, &mut adjoints, output);
            for (slot, &index) in out.iter_mut().zip(wanted) {
                *slot = adjoints.get(index).copied().unwrap_or(T::ZERO);
            }
            Ok(())
        }
    }
}

#[cfg(feature = "alloc")]
pub use growable::GrowTape;
