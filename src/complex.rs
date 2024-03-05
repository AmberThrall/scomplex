use crate::Simplex;
use crate::simplex::Orientation;
use super::errors::OrientationError;
use petgraph::{
    visit::{Bfs, Dfs}, Graph, Undirected
};
use std::collections::HashSet;
use std::ops::{Index,IndexMut};
use std::fmt;
use itertools::*;

/// Pointer to simplex in complex.
pub type SimplexHandle = petgraph::graph::NodeIndex;

/// Struct representing a pair Simplex and its filtration
#[derive(Debug, PartialEq, Clone)]
pub struct FilteredSimplex {
    pub simplex: Simplex,
    pub filtration_value: f32,
}

impl fmt::Display for FilteredSimplex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.simplex, self.filtration_value)
    }
}

/// Struct representing an abstract simplicial complex.
pub struct Complex {
    graph: Graph<FilteredSimplex,f32,Undirected>,
    root: SimplexHandle,
}

impl Complex {
    /// Creates a new abstract simplicial complex
    pub fn new() -> Complex {
        let mut graph = Graph::<FilteredSimplex,f32,Undirected>::new_undirected();
        let root = graph.add_node(FilteredSimplex {
            simplex: Simplex::empty(),
            filtration_value: -1.0,
        });

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
        self.iter().filter(|h| self[*h].simplex.dim() == d).count()
    }

    /// Gets the tree's root (empty simplex)
    pub fn get_root(&self) -> SimplexHandle {
        self.root
    }

    /// Gets a `Vec` of parents of `handle`
    /// Returns an empty `Vec` if `handle is not found.
    pub fn get_parents(&self, handle: SimplexHandle) -> Vec<SimplexHandle> {
        if let Some(sigma) = self.get_simplex(handle) {
            let d = sigma.dim();
            self.graph.neighbors_undirected(handle).filter(|h| self[*h].simplex.dim() == d-1).collect()
        } else {
            Vec::new()
        }
    }

    /// Gets a `Vec` of children of `handle`
    /// Returns an empty `Vec` if `handle is not found.
    pub fn get_children(&self, handle: SimplexHandle) -> Vec<SimplexHandle> {
        if let Some(sigma) = self.get_simplex(handle) {
            let d = sigma.dim();
            self.graph.neighbors_undirected(handle).filter(|h| self[*h].simplex.dim() == d+1).collect()
        } else {
            Vec::new()
        }
    }

    /// Gets a `Vec` of neighbors of `handle`
    /// Returns an empty `Vec` if `handle is not found.
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

    /// Recursively adds a simplex and its faces to the complex, returns a handle to the top simplex.
    /// If the simplex already is part of the complex, it returns the handle to it.
    /// The filtration value defaults to the simplex's dimension.
    pub fn push(&mut self, simplex: Simplex) -> SimplexHandle {
        let d = simplex.dim();
        self.push_filtered(FilteredSimplex {
            simplex,
            filtration_value: d as f32,
        })
    }

    pub fn push_filtered(&mut self, simplex: FilteredSimplex) -> SimplexHandle {
        if let Some(handle) = self.find(&simplex.simplex) {
            return handle;
        }

        if simplex.simplex.dim() > 0 {
            // Add the simplex's faces to the complex
            let mut parents = Vec::new();
            for v in simplex.simplex.vertices() {
                let mut set = simplex.simplex.vertices().clone();
                set.remove(v);

                let h = self.push(Simplex::new(set));
                parents.push(h);
            }

            // Actually add the simplex to the graph
            let handle = self.graph.add_node(simplex);
            for nx in parents {
                self.graph.add_edge(nx, handle, 0.0);
            }
            handle
        } else if simplex.simplex.dim() == 0 {
            let handle = self.graph.add_node(simplex);
            self.graph.add_edge(self.root, handle, 0.0);
            handle
        } else {
            self.root
        }
    }

    /// Looks up a simplex node from a handle
    pub fn get(&self, handle: SimplexHandle) -> Option<&FilteredSimplex> {
        self.graph.node_weight(handle)
    }
    
    pub fn get_mut(&mut self, handle: SimplexHandle) -> Option<&mut FilteredSimplex> {
        self.graph.node_weight_mut(handle)
    }

    /// Looks up a simplex node from a handle, returning just the simplex
    pub fn get_simplex(&self, handle: SimplexHandle) -> Option<&Simplex> {
        self.graph.node_weight(handle).map(|x| &x.simplex)
    }
    
    pub fn get_simplex_mut(&mut self, handle: SimplexHandle) -> Option<&mut Simplex> {
        self.graph.node_weight_mut(handle).map(|x| &mut x.simplex)
    }

    /// Get an iterator over every simplex handle
    pub fn iter(&self) -> petgraph::graph::NodeIndices<petgraph::graph::DefaultIx> {
        self.graph.node_indices()
    }

    /// Looks for a simplex in the complex
    pub fn find(&self, simplex: &Simplex) -> Option<SimplexHandle> {
        self.graph.node_indices().find(|i| &self.graph[*i].simplex == simplex)
    }

    /// Get the complex's dimension, i.e., the largest simplicial dimension in the complex.
    pub fn dim(&self) -> i32 {
        let mut d = -1;
        let mut dfs = Dfs::new(&self.graph, self.root);
        while let Some(nx) = dfs.next(&self.graph) {
            d = d.max(self.get_simplex(nx).unwrap().dim());
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
            let tau = self.get_simplex(nx).unwrap();
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
            let sigma = self.get_simplex_mut(h).unwrap();
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
                let tau = self.get_simplex(face_h).unwrap();
                if tau.orientation() == Orientation::None {
                    continue;
                }

                let tau_orientation = tau.orientation();
                let mut tau_clone = tau.clone();

                let sigma = self.get_simplex_mut(nbr).unwrap();
                sigma.set_orientation(Orientation::Even);
                tau_clone.induced_orientation(&sigma)?;
                if tau_clone.orientation() == tau_orientation {
                    sigma.swap_orientation()?;
                }

                orientations.push(sigma.orientation());
            }

            // Check that we only need one orientation class
            if orientations.len() > 1 {
                for i in 0..orientations.len()-1 {
                    if orientations[i] != orientations[i+1] {
                        self[nbr].simplex.set_orientation(Orientation::None);
                        return Err(OrientationError::Unorientable);
                    }
                }
            }

            // Induce orientation onto unoriented faces and get neighbors
            self._orient_faces(nbr)?;
            for nbr_of_nbr in self.get_neighbors(nbr) {
                if self[nbr_of_nbr].simplex.orientation() == Orientation::None && !neighbors.contains(&nbr_of_nbr) {
                    neighbors.push(nbr_of_nbr);
                }
            }
        }

        Ok(())
    }

    fn _orient_faces(&mut self, root: SimplexHandle) -> Result<(), OrientationError> {
        let sigma = self.get_simplex(root).unwrap().clone();

        for h in self.get_parents(root) {
            let tau = self.get_simplex_mut(h).unwrap();
            if tau.orientation() == Orientation::None {
                tau.induced_orientation(&sigma)?;
                self._orient_faces(h)?;
            }
        }

        Ok(())
    }

    /// Finds all missing `dim` simplices and adds them to the complex.
    /// Sets the filtration value to the max filtration value of its faces.
    pub fn fill_holes(&mut self, dim: i32) {
        if dim <= 0 { return; }

        // Get all (dim+1) combinations of (dim-1) simplices.
        let combinations = self.iter()
            .filter(|h| self[*h].simplex.dim() == dim-1)
            .combinations(dim as usize + 1);

        let mut patches = Vec::new();
        // For each combination, check if the union of simplices is (dim)-dimensional.
        for comb in combinations {
            let mut union: HashSet<usize> = HashSet::new();
            let mut filtration_value: f32 = 0.0;
            for h in comb {
                let sigma = self.get(h).unwrap();
                filtration_value = filtration_value.max(sigma.filtration_value);
                union = union.union(sigma.simplex.vertices())
                    .map(|x| *x)
                    .collect();
            }

            if union.len() == dim as usize + 1 {
                let sigma = Simplex::new(union);
                patches.push(FilteredSimplex {
                    simplex: sigma,
                    filtration_value,
                });
            }
        }

        // Add in the patches.
        for patch in patches {
            self.push_filtered(patch);
        }
    }
}

impl Clone for Complex {
    fn clone(&self) -> Self {
        let graph = self.graph.clone();
        let root = graph.node_indices().find(|i| self.graph[*i].simplex.dim() == -1).unwrap();
        Self {
            graph,
            root,
        }
    }
}

impl Index<SimplexHandle> for Complex {
    type Output = FilteredSimplex;

    fn index(&self, index: SimplexHandle) -> &Self::Output {
        &self.graph[index]
    }
}

impl IndexMut<SimplexHandle> for Complex {
    fn index_mut(&mut self, index: SimplexHandle) -> &mut Self::Output {
        &mut self.graph[index]
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

        assert_eq!(complex[a].simplex.dim(), 0);
        assert_eq!(complex[b].simplex.dim(), 0);
        assert_eq!(complex[ab].simplex.dim(), 1);
        assert_eq!(complex.dim(), 1);
    }

    #[test]
    fn push_recursive() {
        let mut complex = Complex::new();
        let _abc = complex.push(Simplex::from(vec![1, 2, 3]));

        let a = complex.find(&Simplex::from(vec![1])).unwrap();
        let b = complex.find(&Simplex::from(vec![2])).unwrap();
        let c = complex.find(&Simplex::from(vec![3])).unwrap();
        let ab = complex.find(&Simplex::from(vec![1, 2])).unwrap();
        let ac = complex.find(&Simplex::from(vec![1, 3])).unwrap();
        let bc = complex.find(&Simplex::from(vec![2, 3])).unwrap();

        assert_eq!(complex[a].simplex.dim(), 0);
        assert_eq!(complex[b].simplex.dim(), 0);
        assert_eq!(complex[c].simplex.dim(), 0);
        assert_eq!(complex[ab].simplex.dim(), 1);
        assert_eq!(complex[ac].simplex.dim(), 1);
        assert_eq!(complex[bc].simplex.dim(), 1);
        assert_eq!(complex.dim(), 2);
    }

    #[test]
    fn euler_characteristic() {
        let mut complex = Complex::new();
        let _abcd = complex.push(Simplex::from(vec![1, 2, 3, 4]));
        assert_eq!(complex.euler_characteristic(), 4 - 6 + 4 - 1);
    }

    #[test]
    fn skeleton() {
        let mut complex = Complex::new();
        let _abcd = complex.push(Simplex::from(vec![1, 2, 3, 4]));

        let skele2 = complex.skeleton(2);
        let ntri = skele2.num_simplices_by_dimension(2);
        assert_eq!(skele2.dim(), 2);
        assert_eq!(ntri, 4);
    }

    #[test]
    fn neighbors() {
        let mut complex = Complex::new();
        let _abcd = complex.push(Simplex::from(vec![1, 2, 3, 4]));

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
            disk.push(Simplex::from(vec![i,(i + 1) % 5,5]));
        }
        assert!(disk.orient().is_ok());

        let mut mobius_strip = Complex::new();
        mobius_strip.push(Simplex::from(vec![0, 1, 4]));
        mobius_strip.push(Simplex::from(vec![0, 3, 4]));
        mobius_strip.push(Simplex::from(vec![0, 3, 5]));
        mobius_strip.push(Simplex::from(vec![1, 2, 4]));
        mobius_strip.push(Simplex::from(vec![2, 4, 5]));
        mobius_strip.push(Simplex::from(vec![2, 3, 5]));
        assert!(mobius_strip.orient().is_err());
    }
}
