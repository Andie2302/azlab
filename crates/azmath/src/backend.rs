//! Rechen-Backends für die Elementarfunktionen der primitiven Gleitkommatypen.
//!
//! Dieses Modul ist bewusst privat: Der Rest von azmath und alle aufbauenden
//! Crates sehen ausschließlich die Traits
//! [`Elementary`](crate::scalar::Elementary) und
//! [`RealElementary`](crate::scalar::RealElementary).
//! Welches Backend rechnet, entscheidet allein die Feature-Auswahl.
//!
//! ## Neues Backend hinzufügen
//!
//! 1. Optionale Abhängigkeit und Feature in `Cargo.toml` anlegen
//!    (vorher in `todo/03-features-abhaengigkeiten.md` eintragen).
//! 2. Datei `backend/<name>.rs` mit `impl_backend!`-Tabellen für `f32` und
//!    `f64` schreiben.
//! 3. Unten in der Rangfolge eintragen. Es darf immer nur genau ein Backend
//!    übersetzt werden, sonst kollidieren die Trait-Implementierungen.

#[cfg(any(feature = "std", feature = "libm"))]
#[macro_use]
mod table;

// Rangfolge: std > libm.
#[cfg(feature = "std")]
mod std;

#[cfg(all(feature = "libm", not(feature = "std")))]
mod libm;
