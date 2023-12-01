#![no_std]
// Catch documentation errors caused by code changes.
#![deny(rustdoc::broken_intra_doc_links)]

#[cfg(feature = "alloc")]
#[macro_use]
extern crate alloc;

// Re-export ff to make version-matching easier.
pub use ff;

pub mod cofactor;
pub mod curve;
#[cfg(feature = "tests")]
pub mod tests;
pub mod wnaf;

pub use curve::*;
