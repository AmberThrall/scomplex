use super::{
    splx,
    Complex,
    Point,
    PointCloud,
    geom::{Edge, Triangle},
};

/// Complex factory using Delaunay Triangulation using Bowyer-Watson algorithm
/// The resulting simplices represent the indices of the points, i.e., the 1-simplex {0,1} represents the edge `points[0]` -- `points[1]`.
/// Currently only supports a point cloud in R^2.
/// See the `delaunay` example for implementation details.
pub struct Delaunay<'a> {
    points: &'a PointCloud<2>,
}

impl<'a> Delaunay<'a> {
    pub fn new(points: &'a PointCloud<2>) -> Self {
        Delaunay {
            points,
        }
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
            let i0 = self.points.iter().position(|&p| p == tri.v0).unwrap(); 
            let i1 = self.points.iter().position(|&p| p == tri.v1).unwrap(); 
            let i2 = self.points.iter().position(|&p| p == tri.v2).unwrap(); 
            complex.push(splx![i0, i1, i2]);
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
