use crate::Simplex;
use petgraph::{
    Graph,
    Undirected,
    visit::{Bfs, Dfs},
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
        self.handles().filter(|h| self.get(*h).dim() == d).count()
    }

    /// Gets the tree's root (empty simplex)
    pub fn get_root(&self) -> &Simplex {
        self.get(self.root)
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

            if tau.is_face(sigma) && tau.dim() + 1 == d {
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
            for v in &simplex.0 {
                let mut set = simplex.0.clone();
                set.remove(v);
                self.push_recursive(Simplex(set));
            }
        }

        self.push(simplex)
    }

    /// Looks up a simplex node from a handle
    pub fn get(&self, handle: SimplexHandle) -> &Simplex {
        &self.graph[handle]
    }

    /// Get an iterator over every simplex handle
    pub fn handles(&self) -> petgraph::graph::NodeIndices<petgraph::graph::DefaultIx> {
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
        if d < 0 {
            return 0;
        }

        let mut counts = vec![0; d as usize + 1];
        let mut dfs = Dfs::new(&self.graph, self.root);
        while let Some(nx) = dfs.next(&self.graph) {
            let tau_dim = self.get(nx).dim();
            if tau_dim >= 0 {
                counts[tau_dim as usize] += 1;
            }
        }

        let mut x = 0;
        let mut sign = 1i32;
        for i in 0..(d as usize) {
            x += sign * counts[i];
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

    /// Presents the graph structure in DOT format.
    pub fn dot(&self) -> String {
        let dot = petgraph::dot::Dot::with_config(&self.graph, &[petgraph::dot::Config::EdgeNoLabel]);
        format!("{}", dot)
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
    }

    #[test]
    fn euler_characteristic() {
        let mut complex = Complex::new();
        let _abcd = complex.push_recursive(Simplex::from(vec![1, 2, 3, 4]));
        assert_eq!(complex.euler_characteristic(), 2);
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
}
