//! Finite Elemente.
//!
//! - **Elementmatrizen** für die üblichen linearen Elemente: Strecke
//!   ([`line2_stiffness`], …), Dreieck ([`triangle3_stiffness`], …), Viereck
//!   ([`quad4_stiffness`], …), Tetraeder ([`tetra4_stiffness`], …) und das
//!   Dreieck der ebenen Elastizität ([`triangle3_elasticity`]).
//! - [`ElementOperator`] wendet die aus allen Elementen zusammengesetzte
//!   Matrix auf einen Vektor an, ohne sie zu speichern. Zusammen mit den
//!   iterativen Lösern lässt sich so ohne Heap rechnen.
//! - [`assemble_load`] setzt die rechte Seite zusammen,
//!   [`ElementOperator::lift`] trägt vorgegebene Randwerte ein.
//! - Mit dem Feature `alloc` baut `assemble_csr` die Matrix im CSR-Format
//!   auf.
//!
//! Netze kommen vom Aufrufer: eine Liste von Knotenkoordinaten und eine
//! Liste von Elementen, die jeweils ihre Knotennummern nennen. Netzerzeugung
//! gehört nicht zu azmath.

#[cfg(feature = "alloc")]
mod assemble;
mod elements;
mod error;
mod operator;

#[cfg(test)]
mod tests;

#[cfg(feature = "alloc")]
pub use assemble::{CsrOwned, assemble_csr};
pub use elements::{
    PlaneModel, line2_load, line2_mass, line2_stiffness, quad4_mass, quad4_stiffness,
    tetra4_load, tetra4_mass, tetra4_stiffness, triangle3_elasticity, triangle3_load,
    triangle3_mass, triangle3_stiffness,
};
pub use error::FemError;
pub use operator::{ElementOperator, assemble_load};
