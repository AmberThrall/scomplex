use super::{splx, Complex, FilteredSimplex, Point, PointCloud};
use itertools::*;

fn euclidean_distance<const N: usize>(a: &Point<N>, b: &Point<N>) -> f32 {
    let mut s: f32 = 0.0;

    for i in 0..N {
        s += (a[i] - b[i]) * (a[i] - b[i]);
    }

    s.sqrt()
}

/// Complex factory using Vietoris-Rips
/// The resulting simplices represent the indices of the points, i.e., the 1-simplex {0,1} represents the edge `points[0]` -- `points[1]`.
/// See the `rips_complex` example for implementation details.
/// * `max_dim` - determines the the highest dimensional simplices to add to the complex.
/// * `threshold` - only add edges if their distance is less than the threshold
/// * `distance_fn` - distance function used to measure the distance between two points. Defaults to Euclidean distance.
pub struct RipsComplex<'a, const N: usize> {
    points: &'a PointCloud<N>,
    threshold: f32,
    max_dim: i32,
    distance_fn: Box<dyn Fn(&Point<N>,&Point<N>) -> f32>,
}

impl<'a, const N: usize> RipsComplex<'a,N> {
    pub fn new(points: &'a PointCloud<N>, threshold: f32) -> Self {
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

    pub fn distance_fn(mut self, f: impl Fn(&Point<N>,&Point<N>) -> f32 + 'static) -> Self {
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
                    complex.push_filtered_simplex(FilteredSimplex {
                        simplex: splx![comb[0], comb[1]],
                        filtration_value: dist,
                    });
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
