use crate::Simplex;
use crate::simplex::Orientation;
use super::errors::OrientationError;
use petgraph::{
    visit::{Bfs, Dfs}, Graph, Undirected
};

/// Pointer to simplex in complex.
pub type SimplexHandle = petgraph::graph::NodeIndex;

/// Struct representing an abstract simplicial complex.
pub struct Complex {
    graph: Graph<Simplex,f32,Undirected>,
    root: SimplexHandle,
}

impl Complex {
    /// Creates a new abstract simplicial complex
    pub fn new() -> Complex {
        let mut graph = Graph::<Simplex,f32,Undirected>::new_undirected();
        let root = graph.add_node(Simplex::empty());

        Complex {
            graph,
            root,
        }
    }

    /// Returns the number of simplices
    pub fn num_simplices(&self) -> usize {
        self.graph.node_count()
    }

    /// Returns the number of simplicies by dimension
    pub fn num_simplices_by_dimension(&self, d: i32) -> usize {
        self.iter().filter(|h| self.get(*h).dim() == d).count()
    }

    /// Gets the tree's root (empty simplex)
    pub fn get_root(&self) -> &Simplex {
        self.get(self.root)
    }

    pub fn get_root_mut(&mut self) -> &mut Simplex {
        self.get_mut(self.root)
    }

    /// Gets a `Vec` of parents of `handle`
    pub fn get_parents(&self, handle: SimplexHandle) -> Vec<SimplexHandle> {
        let d = self.get(handle).dim();
        self.graph.neighbors_undirected(handle).filter(|h| self.get(*h).dim() == d-1).collect()
    }

    /// Gets a `Vec` of children of `handle`
    pub fn get_children(&self, handle: SimplexHandle) -> Vec<SimplexHandle> {
        let d = self.get(handle).dim();
        self.graph.neighbors_undirected(handle).filter(|h| self.get(*h).dim() == d+1).collect()
    }

    /// Gets a `Vec` of neighbors of `handle`
    pub fn get_neighbors(&self, handle: SimplexHandle) -> Vec<SimplexHandle> {
        let parents = self.get_parents(handle);

        let mut neighbors = Vec::new();
        for parent in parents {
            for child in self.get_children(parent) {
                if child != handle {
                    neighbors.push(child);
                }
            }
        }

        neighbors
    }

    /// Append a simplex to the complex, returning a handle to the simplex.
    /// If the simplex already is part of the complex, it returns the handle to it.
    pub fn push(&mut self, simplex: Simplex) -> SimplexHandle {
        if let Some(handle) = self.find(&simplex) {
            return handle;
        }

        // Actually add the simplex to the graph
        let handle = self.graph.add_node(simplex);
        let sigma = self.get(handle);
        let d = sigma.dim();

        let mut bfs = Bfs::new(&self.graph, self.root);
        let mut parents = Vec::new();

        while let Some(nx) = bfs.next(&self.graph) {
            let tau = self.get(nx);
            if tau.dim() + 1 > d {
                break;
            } else if tau.is_face(sigma) && tau.dim() + 1 == d {
                parents.push(nx);
            }
        }
        
        for nx in parents {
            self.graph.add_edge(nx, handle, 0.0);
        }

        handle
    }

    /// Recursively adds a simplex and its faces to the complex, returns a handle to the top simplex.
    pub fn push_recursive(&mut self, simplex: Simplex) -> SimplexHandle {
        // Add the simplex's faces
        if simplex.dim() > 0 {
            for v in simplex.vertices() {
                let mut set = simplex.vertices().clone();
                set.remove(v);
                self.push_recursive(Simplex::new(set));
            }
        }

        self.push(simplex)
    }

    /// Looks up a simplex node from a handle
    pub fn get(&self, handle: SimplexHandle) -> &Simplex {
        &self.graph[handle]
    }
    
    pub fn get_mut(&mut self, handle: SimplexHandle) -> &mut Simplex {
        &mut self.graph[handle]
    }

    /// Get an iterator over every simplex handle
    pub fn iter(&self) -> petgraph::graph::NodeIndices<petgraph::graph::DefaultIx> {
        self.graph.node_indices()
    }

    /// Looks for a simplex in the complex
    pub fn find(&self, simplex: &Simplex) -> Option<SimplexHandle> {
        self.graph.node_indices().find(|i| &self.graph[*i] == simplex)
    }

    /// Get the complex's dimension, i.e., the largest simplicial dimension in the complex.
    pub fn dim(&self) -> i32 {
        let mut d = -1;
        let mut dfs = Dfs::new(&self.graph, self.root);
        while let Some(nx) = dfs.next(&self.graph) {
            d = d.max(self.get(nx).dim());
        }

        d
    }

    /// Calculate the Euler characteristic
    pub fn euler_characteristic(&self) -> i32 {
        let d = self.dim();
        if d <= 0 {
            return 0;
        }

        let mut x = 0;
        let mut sign = 1i32;
        for i in 0..(d as usize)+1 {
            x += sign * (self.num_simplices_by_dimension(i as i32) as i32);
            sign *= -1;
        }
        x
    }

    /// Get the j-skeleton of the complex, i.e., the subcomplex of simplices with dim <= j
    pub fn skeleton(&self, j: i32) -> Complex {
        let mut subcomplex = Complex::new();

        let mut bfs = Bfs::new(&self.graph, self.root);
        while let Some(nx) = bfs.next(&self.graph) {
            let tau = self.get(nx);
            if tau.dim() > j {
                break;
            }

            subcomplex.push(tau.clone());
        }

        subcomplex
    }

    /// Creates a Hasse diagram in DOT format.
    pub fn hasse(&self) -> String {
        let dot = petgraph::dot::Dot::with_config(&self.graph, &[petgraph::dot::Config::EdgeNoLabel]);
        format!("{}", dot)
    }

    /// Attempts to orient the simplices.
    pub fn orient(&mut self) -> Result<(), OrientationError> {
        // 1. Set the orientation of a random d-simplex and induce its orientation onto its faces, then
        //    get its neighbors
        let d = self.dim();
        let mut neighbors = Vec::new();
        for h in self.iter() {
            let sigma = self.get_mut(h);
            sigma.set_orientation(Orientation::None);

            if sigma.dim() == d && neighbors.is_empty() {
                sigma.set_orientation(Orientation::Even);
                self._orient_faces(h)?;
                neighbors = self.get_neighbors(h);
            }
        }

        // 2. For each neighbor:
        //   a. Set the orientation such that the induced orientation on its vfaces is opposite
        //      their current orientation
        //   b. Induce orientation onto non-oriented faces
        //   c. Add each unoriented neighbor to the neighbors list
        while let Some(nbr) = neighbors.pop() {
            let mut orientations: Vec<Orientation> = Vec::new();

            // Determine the needed orientations
            for face_h in self.get_parents(nbr) {
                let tau = self.get(face_h);
                if tau.orientation() == Orientation::None {
                    continue;
                }

                let tau_orientation = tau.orientation();
                let mut tau_clone = tau.clone();

                let sigma = self.get_mut(nbr);
                sigma.set_orientation(Orientation::Even);
                tau_clone.induced_orientation(sigma)?;
                if tau_clone.orientation() == tau_orientation {
                    sigma.swap_orientation()?;
                }

                orientations.push(sigma.orientation());
            }

            // Check that we only need one orientation class
            if orientations.len() > 1 {
                for i in 0..orientations.len()-1 {
                    if orientations[i] != orientations[i+1] {
                        self.get_mut(nbr).set_orientation(Orientation::None);
                        return Err(OrientationError::Unorientable);
                    }
                }
            }

            // Induce orientation onto unoriented faces and get neighbors
            self._orient_faces(nbr)?;
            for nbr_of_nbr in self.get_neighbors(nbr) {
                if self.get(nbr_of_nbr).orientation() == Orientation::None && !neighbors.contains(&nbr_of_nbr) {
                    neighbors.push(nbr_of_nbr);
                }
            }
        }

        Ok(())
    }

    fn _orient_faces(&mut self, root: SimplexHandle) -> Result<(), OrientationError> {
        let sigma = self.get(root).clone();

        for h in self.get_parents(root) {
            let tau = self.get_mut(h);
            if tau.orientation() == Orientation::None {
                tau.induced_orientation(&sigma)?;
                self._orient_faces(h)?;
            }
        }

        Ok(())
    }
}

impl Clone for Complex {
    fn clone(&self) -> Self {
        let graph = self.graph.clone();
        let root = graph.node_indices().find(|i| self.graph[*i].dim() == -1).unwrap();
        Self {
            graph,
            root,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn push_get() {
        let mut complex = Complex::new();
        let a = complex.push(Simplex::from(vec![1]));
        let b = complex.push(Simplex::from(vec![2]));
        let ab = complex.push(Simplex::from(vec![1, 2]));

        assert_eq!(complex.get(a).dim(), 0);
        assert_eq!(complex.get(b).dim(), 0);
        assert_eq!(complex.get(ab).dim(), 1);
        assert_eq!(complex.dim(), 1);
    }

    #[test]
    fn push_recursive() {
        let mut complex = Complex::new();
        let _abc = complex.push_recursive(Simplex::from(vec![1, 2, 3]));

        let a = complex.find(&Simplex::from(vec![1])).unwrap();
        let b = complex.find(&Simplex::from(vec![2])).unwrap();
        let c = complex.find(&Simplex::from(vec![3])).unwrap();
        let ab = complex.find(&Simplex::from(vec![1, 2])).unwrap();
        let ac = complex.find(&Simplex::from(vec![1, 3])).unwrap();
        let bc = complex.find(&Simplex::from(vec![2, 3])).unwrap();

        assert_eq!(complex.get(a).dim(), 0);
        assert_eq!(complex.get(b).dim(), 0);
        assert_eq!(complex.get(c).dim(), 0);
        assert_eq!(complex.get(ab).dim(), 1);
        assert_eq!(complex.get(ac).dim(), 1);
        assert_eq!(complex.get(bc).dim(), 1);
        assert_eq!(complex.dim(), 2);
        // println!("{}", complex.hasse());
    }

    #[test]
    fn euler_characteristic() {
        let mut complex = Complex::new();
        let _abcd = complex.push_recursive(Simplex::from(vec![1, 2, 3, 4]));
        assert_eq!(complex.euler_characteristic(), 4 - 6 + 4 - 1);
    }

    #[test]
    fn skeleton() {
        let mut complex = Complex::new();
        let _abcd = complex.push_recursive(Simplex::from(vec![1, 2, 3, 4]));

        let skele2 = complex.skeleton(2);
        let ntri = skele2.num_simplices_by_dimension(2);
        assert_eq!(skele2.dim(), 2);
        assert_eq!(ntri, 4);
    }

    #[test]
    fn neighbors() {
        let mut complex = Complex::new();
        let _abcd = complex.push_recursive(Simplex::from(vec![1, 2, 3, 4]));

        let abc = complex.find(&Simplex::from(vec![1,2,3])).unwrap();
        let abd = complex.find(&Simplex::from(vec![1,2,4])).unwrap();
        let acd = complex.find(&Simplex::from(vec![1,3,4])).unwrap();
        let bcd = complex.find(&Simplex::from(vec![2,3,4])).unwrap();

        let neighbors = complex.get_neighbors(abc);
        assert!(neighbors.contains(&abd));
        assert!(neighbors.contains(&acd));
        assert!(neighbors.contains(&bcd));
    }

    #[test]
    fn orient() {
        let mut disk = Complex::new();
        for i in 0..5 {
            disk.push_recursive(Simplex::from(vec![i,(i + 1) % 5,5]));
        }
        assert!(disk.orient().is_ok());

        let mut mobius_strip = Complex::new();
        mobius_strip.push_recursive(Simplex::from(vec![0, 1, 4]));
        mobius_strip.push_recursive(Simplex::from(vec![0, 3, 4]));
        mobius_strip.push_recursive(Simplex::from(vec![0, 3, 5]));
        mobius_strip.push_recursive(Simplex::from(vec![1, 2, 4]));
        mobius_strip.push_recursive(Simplex::from(vec![2, 4, 5]));
        mobius_strip.push_recursive(Simplex::from(vec![2, 3, 5]));
        assert!(mobius_strip.orient().is_err());
    }
}
