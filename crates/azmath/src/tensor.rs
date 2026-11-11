//! Tensoren beliebiger Stufe.
//!
//! Die Stufe (Anzahl der Indizes) steht im Typ, die Form erst zur Laufzeit.
//! [`TensorView`] und [`TensorViewMut`] sind Sichten auf Puffer des
//! Aufrufers und brauchen keinen Heap. Gespeichert wird so, dass der letzte
//! Index am schnellsten läuft.
//!
//! Operationen, die die Stufe ändern ([`outer`], [`dot`], [`contract`]),
//! schreiben ihr Ergebnis in einen Puffer des Aufrufers; daraus macht man
//! mit [`TensorView::new`] wieder einen Tensor der passenden Stufe.
//!
//! Mit dem Feature `alloc` gibt es zusätzlich den besitzenden `DTensor`.

mod ops;
#[cfg(feature = "alloc")]
mod owned;
mod view;

#[cfg(test)]
mod tests;

pub use ops::{contract, dot, outer};
#[cfg(feature = "alloc")]
pub use owned::DTensor;
pub use view::{TensorView, TensorViewMut};
