use super::{
    splx,
    FilteredSimplex,
    Complex,
    Point,
    PointCloud,
    geom::{Edge, Triangle},
};

/// Complex factory using Alpha Complex via Bowyer-Watson algorithm
/// The resulting simplices represent the indices of the points, i.e., the 1-simplex {0,1} represents the edge `points[0]` -- `points[1]`.
/// Sets the filtration values for 1-simplices to the edge length, and the filtration values for
/// 2-simplices as the max of filtration values among faces. Only adds simplices whose filtration
/// values is <= `alpha`.
/// Currently only supports a point cloud in R^2.
/// See the `delaunay_complex` example for implementation details.
pub struct AlphaComplex<'a> {
    points: &'a PointCloud<2>,
    alpha: f32,
}

impl<'a> AlphaComplex<'a> {
    pub fn delaunay(points: &'a PointCloud<2>) -> Self {
        AlphaComplex::new(points, f32::MAX)
    }
 
    pub fn new(points: &'a PointCloud<2>, alpha: f32) -> Self {
        AlphaComplex {
            points,
            alpha,
        }
    }

    pub fn alpha(mut self, alpha: f32) -> Self {
        self.alpha = alpha;
        self
    }

    pub fn build(self) -> Complex {
        // Create triangle bounding every point
        let st = self.super_triangle();
        let mut triangles = vec![st.clone()];

        // Triangulate each point 
        for pt in self.points {
            triangles = self.add_point(&pt, triangles);
        }

        // Remove triangles that share an edge with the super triangle.
        let triangles = triangles.iter().filter(|tri| {
            !(tri.v0 == st.v0 || tri.v0 == st.v1 || tri.v0 == st.v2 ||
              tri.v1 == st.v0 || tri.v1 == st.v1 || tri.v1 == st.v2 ||
              tri.v2 == st.v0 || tri.v2 == st.v1 || tri.v2 == st.v2)
        });
        
        // Construct the complex.
        let mut complex = Complex::new();
        for tri in triangles {
            let v0tov1 = crate::geom::euclidean_distance(&tri.v0, &tri.v1) / 2.0;
            let v0tov2 = crate::geom::euclidean_distance(&tri.v0, &tri.v2) / 2.0;
            let v1tov2 = crate::geom::euclidean_distance(&tri.v1, &tri.v2) / 2.0;

            let max_dist = v0tov1.max(v0tov2).max(v1tov2);

            let i0 = self.points.iter().position(|&p| p == tri.v0).unwrap(); 
            let i1 = self.points.iter().position(|&p| p == tri.v1).unwrap(); 
            let i2 = self.points.iter().position(|&p| p == tri.v2).unwrap(); 

            if v0tov1 <= self.alpha {
                complex.push_filtered_simplex(FilteredSimplex {
                    simplex: splx![i0, i1],
                    filtration_value: v0tov1,
                });
            }

            if v0tov2 <= self.alpha {
                complex.push_filtered_simplex(FilteredSimplex {
                    simplex: splx![i0, i2],
                    filtration_value: v0tov2,
                });
            }

            if v1tov2 <= self.alpha {
                complex.push_filtered_simplex(FilteredSimplex {
                    simplex: splx![i1, i2],
                    filtration_value: v1tov2,
                });
            }

            if max_dist <= self.alpha {
                complex.push_filtered_simplex(FilteredSimplex {
                    simplex: splx![i0, i1, i2],
                    filtration_value: max_dist,
                });
            }

        }
        complex
    }

    fn super_triangle(&self) -> Triangle {
        // Find the maximum and minimum coordinates and find a triangle that encloses all points.
        let (mut min_x, mut min_y) = (f32::MAX, f32::MAX);
        let (mut max_x, mut max_y) = (f32::MIN, f32::MIN);

        for pt in self.points {
            min_x = min_x.min(pt[0]);
            min_y = min_y.min(pt[1]);
            max_x = max_x.max(pt[0]);
            max_y = max_y.max(pt[1]);
        }

        let dx = (max_x - min_x) * 10.0;
        let dy = (max_y - min_y) * 10.0;

        let v0 = [min_x - dx, min_y - dy * 3.0];
        let v1 = [min_x - dx, max_y + dy];
        let v2 = [max_x + dx * 3.0, max_y + dy];

        Triangle::new(v0, v1, v2)
    }

    fn add_point(&self, pt: &Point<2>, mut triangles: Vec<Triangle>) -> Vec<Triangle> {
        let mut edges = Vec::new();

        // Remove any triangles whose circumcricle contain the vertex
        triangles.retain(|tri| {
            if tri.circumcircle.contains(pt) {
                edges.push(Edge::new(tri.v0, tri.v1));
                edges.push(Edge::new(tri.v1, tri.v2));
                edges.push(Edge::new(tri.v2, tri.v0));
                false
            } else {
                true 
           }
        });

        // Get the unique edges
        let mut unique_edges = Vec::new();
        for i in 0..edges.len() {
            let count = edges.iter().filter(|e| e == &&edges[i]).count();

            if count == 1 {
                unique_edges.push(edges[i]);
            }
        }

        // Add a triangle for each unique edge
        for edge in unique_edges {
            let tri = Triangle::new(edge.v0, edge.v1, pt.clone()); 
            triangles.push(tri);
        }

        triangles
    }
}
