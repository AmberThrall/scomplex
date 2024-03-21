#![allow(dead_code)]
pub mod errors;
pub mod simplex;
pub mod complex;
pub mod rips_complex;
pub mod alpha_complex;
pub mod geom;
pub mod chain;

pub use errors::*;
pub use simplex::Simplex;
pub use complex::Complex;
pub use complex::FilteredSimplex;
pub use rips_complex::RipsComplex;
pub use alpha_complex::AlphaComplex;
pub use chain::Chain;

pub type Point<const N: usize> = [f32; N];
pub type PointCloud<const N: usize> = Vec<Point<N>>;

