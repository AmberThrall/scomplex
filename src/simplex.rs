use std::collections::{HashMap, HashSet};
use std::fmt;
use super::errors::OrientationError;

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

/// Possible Orientations
#[derive(Debug,PartialEq,Eq)]
pub enum Orientation {
    None,
    Even,
    Odd,
}

/// Struct representing an abstract simplex.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Simplex {
    vertices: HashSet<usize>,
    orientation: Vec<usize>,
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
            orientation: Vec::new(),
        }
    }

    /// Returns the vertices `HashSet`
    pub fn vertices(&self) -> &HashSet<usize> {
        &self.vertices
    }

    /// Gets the current orientation of the simplex.
    pub fn orientation(&self) -> Orientation {
        if self.orientation.is_empty() {
            return Orientation::None;
        }
        if self.dim() <= 0 {
            return Orientation::Even;
        }

        // Source: https://www.geeksforgeeks.org/minimum-number-swaps-required-sort-array/
        let mut nums = self.orientation.clone();
        let mut map = HashMap::new();
        for i in 0..nums.len() {
            map.insert(nums[i], i);
        }
        nums.sort();

        let mut visited = vec![false; nums.len()];
        let mut num_swaps = 0;
        for i in 0..nums.len() {
            if visited[i] || map[&nums[i]] == i {
                continue;
            }

            let mut j = i;
            let mut cycle_size = 0;
            while visited[j] == false {
                visited[j] = true;
                j = map[&nums[j]];
                cycle_size += 1;
            }

            if cycle_size > 0 {
                num_swaps += cycle_size - 1;
            }
        }


        if num_swaps % 2 == 0 {
            Orientation::Even
        } else {
            Orientation::Odd
        }
    }

    /// Sets the orientation of the simplex.
    pub fn set_orientation(&mut self, orientation: Orientation) {
        if orientation == Orientation::None {
            self.orientation = Vec::new();
            return;
        }

        self.orientation = self.vertices.iter()
            .map(|x| *x)
            .collect();
        self.orientation.sort();
        if orientation == Orientation::Odd {
            self.swap_orientation().unwrap();
        }
    }

    /// Swaps the orientation of the simplex. Returns an error is the simplex is not oriented.
    pub fn swap_orientation(&mut self) -> Result<(),OrientationError> {
        if self.orientation.is_empty() {
            Err(OrientationError::CannotSwapOrientation)
        } else if self.dim() <= 0 {
            Ok(())
        } else {
            let tmp = self.orientation[0];
            self.orientation[0] = self.orientation[1];
            self.orientation[1] = tmp;
            Ok(())
        }
    }

    /// Gets the induced orientation of the simplex from its coface `other`
    pub fn induced_orientation(&mut self, other: &Simplex) -> Result<(), OrientationError> {
        if !self.is_face(other) {
            return Err(OrientationError::NotACoface);
        }
        if self.dim() < 0 {
            return Ok(());
        }

        self.orientation = Vec::with_capacity(self.dim() as usize + 1);
        let mut odd_idx = false;
        for (i, v) in other.orientation.iter().enumerate() {
            if !self.vertices.contains(v) {
                odd_idx = i % 2 == 1;
            } else {
                self.orientation.push(*v);
            }
        }

        if odd_idx{
            self.swap_orientation().unwrap();
        }
        Ok(())
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
        if self.orientation.is_empty() {
            write!(f, "{:?}", self.vertices)
        } else {
            write!(f, "{:?}", self.orientation)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::simplex::Orientation;

    #[test]
    fn simplex() {
        let simplex1 = splx![1,2,3];
        let simplex2 = splx![3,2,1];
        let face = splx![1,2];
        assert_eq!(simplex1, simplex2);
        assert!(face.is_face(&simplex1));
        assert!(simplex1.is_coface(&face));
    }

    #[test]
    fn orientation() {
        let mut sigma = splx![1,2,3];

        assert_eq!(sigma.orientation(), Orientation::None);

        // sigma = [1,2,3]
        sigma.set_orientation(Orientation::Even);
        assert_eq!(sigma.orientation(), Orientation::Even);

        // sigma = [2,1,3]
        assert!(sigma.swap_orientation().is_ok());
        assert_eq!(sigma.orientation(), Orientation::Odd);

        // tau = [3,2]
        let mut tau = splx![2,3];
        assert!(tau.induced_orientation(&sigma).is_ok());
        assert_eq!(sigma.orientation(), Orientation::Odd);
    }
}
