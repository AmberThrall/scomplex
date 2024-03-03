use super::{splx, Complex, Point, PointCloud};
use itertools::*;

fn euclidean_distance(a: &Point, b: &Point) -> f32 {
    let mut s: f32 = 0.0;

    for i in 0..a.len() {
        s += (a[i] - b[i]) * (a[i] - b[i]);
    }

    s.sqrt()
}

pub struct RipsComplex<'a> {
    points: &'a PointCloud,
    threshold: f32,
    max_dim: i32,
    distance_fn: Box<dyn Fn(&Point,&Point) -> f32>,
}

impl<'a> RipsComplex<'a> {
    pub fn new(points: &'a PointCloud, threshold: f32) -> Self {
        RipsComplex {
            points,
            threshold,
            max_dim: 1,
            distance_fn: Box::new(euclidean_distance),
        }
    }

    pub fn max_dim(mut self, dim: i32) -> Self {
        self.max_dim = dim;
        self
    }

    pub fn threshold(mut self, threshold: f32) -> Self {
        self.threshold = threshold;
        self
    }

    pub fn distance_fn(mut self, f: impl Fn(&Point,&Point) -> f32 + 'static) -> Self {
        self.distance_fn = Box::new(f);
        self
    }

    pub fn build(self) -> Complex {
        let mut complex = Complex::new();

        // Add each vertex as a 0-simplex
        for i in 0..self.points.len() {
            complex.push(splx![i]);
        }

        // Add every edge ab if the distance from a to b is less than radius
        if self.max_dim > 0 {
            for comb in (0..self.points.len()).combinations(2) {
                let dist = (self.distance_fn)(&self.points[comb[0]], &self.points[comb[1]]);
                if dist < self.threshold {
                    complex.push(splx![comb[0], comb[1]]);
                }
            }
        }

        if self.max_dim > 1 {
            for k in 2..self.max_dim+1 {
                complex.fill_holes(k);
            }
        }

        complex
    }
}