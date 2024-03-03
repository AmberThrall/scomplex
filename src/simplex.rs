use std::collections::HashSet;
use std::fmt;

/// Creates a `Simplex` with vertices given by the arguemnts
/// 
/// ```
/// use scomplex::splx;
/// let simplex = splx![1,2,3];
/// assert_eq!(simplex.dim(), 2);
/// ```
#[macro_export]
macro_rules! splx {
    ( $( $x:expr ),* ) => {
        {
            let mut temp_set = std::collections::HashSet::new();
            $(
                temp_set.insert($x);
            )*
            $crate::simplex::Simplex::new(temp_set)
        }
    }
}

/// Struct representing an abstract simplex.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Simplex {
    pub vertices: HashSet<usize>,
}

impl Simplex {
    /// Creates an empty simplex.
    pub fn empty() -> Simplex {
        Simplex::new(HashSet::new())
    }

    /// Creates an unoriented simplex.
    pub fn new(vertices: HashSet<usize>) -> Simplex {
        Simplex {
            vertices,
        }
    }

    /// Returns the dimension of the simplex, i.e., the number of vertices minus 1.
    /// An empty simplex is treated as having dimension -1.
    pub fn dim(&self) -> i32 {
        (self.vertices.len() as i32) - 1
    }

    /// Returns `true` if the simplex is a face of another simplex `other`.
    pub fn is_face(&self, other: &Simplex) -> bool {
        self.vertices.is_subset(&other.vertices)
    }

    /// Returns `true` if the simplex is a coface of another simplex `other`.
    pub fn is_coface(&self, other: &Simplex) -> bool {
        other.is_face(&self)
    }
}

impl From<Vec<usize>> for Simplex {
    fn from(value: Vec<usize>) -> Simplex {
        let mut set = HashSet::new();
        for v in value {
            set.insert(v);
        }
        Simplex::new(set)
    }
}

impl fmt::Display for Simplex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.vertices)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn simplex() {
        let simplex1 = splx![1,2,3];
        let simplex2 = splx![3,2,1];
        let face = splx![1,2];
        assert_eq!(simplex1, simplex2);
        assert!(face.is_face(&simplex1));
        assert!(simplex1.is_coface(&face));
    }
}
