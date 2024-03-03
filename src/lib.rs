#![allow(dead_code)]
pub mod errors;
pub mod simplex;
pub mod complex;
pub mod rips_complex;

pub use errors::*;
pub use simplex::Simplex;
pub use complex::Complex;
pub use rips_complex::RipsComplex;

pub type Point = Vec<f32>;
pub type PointCloud = Vec<Point>;