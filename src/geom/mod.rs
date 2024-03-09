extern crate nalgebra as na;
use super::Point;

pub fn squared_euclidean_distance<const N: usize>(a: &Point<N>, b: &Point<N>) -> f32 {
    let mut s: f32 = 0.0;

    for i in 0..N {
        s += (a[i] - b[i]) * (a[i] - b[i]);
    }

    s
}

pub fn euclidean_distance<const N: usize>(a: &Point<N>, b: &Point<N>) -> f32 {
    squared_euclidean_distance(a, b).sqrt()
}

#[derive(Debug,PartialEq,Clone,Copy)]
pub struct Circle {
    pub center: Point<2>,
    pub radius: f32,
}

impl Circle {
    pub fn new(center: Point<2>, radius: f32) -> Self {
        Circle {
            center, radius
        }
    }

    pub fn contains(&self, pt: &Point<2>) -> bool {
        let dx = self.center[0] - pt[0];
        let dy = self.center[1] - pt[1];
        let r = dx * dx + dy * dy;
        r <= self.radius * self.radius
    }
}

#[derive(Debug,Clone,Copy)]
pub struct Edge {
    pub v0: Point<2>,
    pub v1: Point<2>,
}

impl Edge {
    pub fn new(v0: Point<2>, v1: Point<2>) -> Self {
        Edge {
            v0, v1 
        }
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.v0 == other.v0 && self.v1 == other.v1) || 
            (self.v0 == other.v1 && self.v1 == other.v0)
    }
}

#[derive(Debug,Clone,Copy)]
pub struct Triangle {
    pub v0: Point<2>,
    pub v1: Point<2>,
    pub v2: Point<2>,
    pub circumcircle: Circle,
}

impl Triangle {
    pub fn new(v0: Point<2>, v1: Point<2>, v2: Point<2>) -> Self {
       let mut tri = Triangle {
            v0, v1, v2,
            circumcircle: Circle {
                center: v0,
                radius: 0.0
            },
       };

       tri.compute_circumcircle();
       tri
    }

    pub fn compute_circumcircle(&mut self) {
        // See https://ics.uci.edu/~eppstein/junkyard/circumcenter.html
        // Original SoE:
        //   (x-a)^2 + (y-b)^2 = r^2
        // which expands out to be
        //   x^2 - 2ax + a^2 + y^2 - 2by + b^2 = r^2
        // let q = r^2 - a^2 - b^2 and rearrange to get 
        //   2xa + 2yb + q = x^2 + y^2
        // This is a linear system of 3 unknowns (a, b, & q) with 3 equations (3 points on circle)
        // which can be solved via Cramer's rule
        
        let m = na::Matrix3::new(2.0 * self.v0[0], 2.0 * self.v0[1], 1.0,
                                 2.0 * self.v1[0], 2.0 * self.v1[1], 1.0,
                                 2.0 * self.v2[0], 2.0 * self.v2[1], 1.0);

        let y = na::Vector3::new(self.v0[0] * self.v0[0] + self.v0[1] * self.v0[1],
                                 self.v1[0] * self.v1[0] + self.v1[1] * self.v1[1],
                                 self.v2[0] * self.v2[0] + self.v2[1] * self.v2[1]);

        let m1 = na::Matrix3::from_columns(&[y, m.column(1).into(), m.column(2).into()]);
        let m2 = na::Matrix3::from_columns(&[m.column(0).into(), y, m.column(2).into()]);

        let denom = m.determinant();
        let a = m1.determinant() / denom;
        let b = m2.determinant() / denom;
        let r = ((self.v0[0] - a) * (self.v0[0] - a) + (self.v0[1] - b) * (self.v0[1] - b)).sqrt();

        self.circumcircle = Circle {
            center: [a,b],
            radius: r,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circumcircle() {
        let tri = Triangle::new([3.0,2.0], [1.0,-1.0], [-3.0,1.0]);
        let center = [-0.0625, 1.875];
        let radius = euclidean_distance(&[3.0,2.0], &center);

        assert_eq!(tri.circumcircle.center, center);
        assert!((tri.circumcircle.radius - radius).abs() < 1e-13f32);
    }
} 
