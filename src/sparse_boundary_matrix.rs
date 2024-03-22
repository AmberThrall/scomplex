use crate::{
    Complex,
    complex::SimplexHandle
};
use std::collections::{HashSet, HashMap};
use std::fmt;

/// Represents a boundary matrix of Z_2
/// only stores shape and non-zero entries.
#[derive(Debug, PartialEq)]
pub struct SparseBoundaryMatrix {
    shape: (usize, usize),
    pub entries: HashSet<(usize, usize)>,
}

impl SparseBoundaryMatrix {
    /// Construct the p-th boundary matrix
    pub fn new(p: i32, complex: &Complex) -> SparseBoundaryMatrix {
        // associate a column for each p-simplex and a row for each (p-1)-simplex 
        let mut col_lookup: HashMap<SimplexHandle, usize> = HashMap::new(); 
        let mut row_lookup: HashMap<SimplexHandle, usize> = HashMap::new(); 

        for h in complex.iter() {
            let dim = complex.get_simplex(h).unwrap().dim();
            if dim == p {
                col_lookup.insert(h, col_lookup.len());
            } else if dim == p-1 {
                row_lookup.insert(h, row_lookup.len());
            }
        }

        // add the entries
        let shape = (row_lookup.len(), col_lookup.len());
        let mut entries = HashSet::new();
        for (h,j) in col_lookup.iter() {
            for face in complex.get_parents(*h) {
                let i = row_lookup.get(&face).unwrap();
                entries.insert((*i,*j));
            }
        }

        SparseBoundaryMatrix { 
            shape,
            entries 
        }
    }

    /// Returns the matrix shape (#rows, #cols)
    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    /// Gets an entry from the matrix
    pub fn get(&self, i: usize, j: usize) -> u32 {
        match self.entries.contains(&(i,j)) {
            true => 1,
            false => 0,
        }
    }

    /// Sets the matrix entry at (i,j) to v % 2
    pub fn set(&mut self, i: usize, j: usize, v: u32) {
        if v % 2 == 1 {
            self.entries.insert((i,j));
        } else {
            self.entries.remove(&(i,j));
        }
    }

    /// Swaps two rows: Ri <-> Rj
    pub fn swap_rows(&mut self, i: usize, j: usize) {
        // Get copies of the rows
        let rowi: Vec<(usize,usize)> = self.entries.iter()
            .filter(|(r,_c)| *r == i)
            .map(|(r,c)| (*r,*c))
            .collect();
        let rowj: Vec<(usize,usize)> = self.entries.iter()
            .filter(|(r,_c)| *r == j)
            .map(|(r,c)| (*r,*c))
            .collect();

        // Clear both rows
        for (r,c) in rowi.iter() { self.set(*r,*c,0); }
        for (r,c) in rowj.iter() { self.set(*r,*c,0); }

        // Fill in rows
        for (_r,c) in rowi { self.set(j, c, 1); }
        for (_r,c) in rowj { self.set(i, c, 1); }
    }

    /// Swaps two columns: Ri <-> Rj
    pub fn swap_columns(&mut self, i: usize, j: usize) {
        // Get copies of the columns
        let coli: Vec<(usize,usize)> = self.entries.iter()
            .filter(|(_r,c)| *c == i)
            .map(|(r,c)| (*r,*c))
            .collect();
        let colj: Vec<(usize,usize)> = self.entries.iter()
            .filter(|(_r,c)| *c == j)
            .map(|(r,c)| (*r,*c))
            .collect();

        // Clear both rows
        for (r,c) in coli.iter() { self.set(*r,*c,0); }
        for (r,c) in colj.iter() { self.set(*r,*c,0); }

        // Fill in rows
        for (r,_c) in coli { self.set(r, j, 1); }
        for (r,_c) in colj { self.set(r, i, 1); }
    }
        
    /// Performs row replacement: R_i <- R_i + R_j
    pub fn add_row(&mut self, i: usize, j: usize) { 
        let rowj: Vec<(usize,usize)> = self.entries.iter()
            .filter(|(r,_c)| *r == j)
            .map(|(r,c)| (*r,*c))
            .collect();

        for (_r,c) in rowj {
            let v = self.get(i,c) + 1;
            self.set(i, c, v); 
        }
    }
    
    /// Performs column replacement: C_i <- C_i + C_j
    pub fn add_column(&mut self, i: usize, j: usize) { 
        let colj: Vec<(usize,usize)> = self.entries.iter()
            .filter(|(_r,c)| *c == j)
            .map(|(r,c)| (*r,*c))
            .collect();

        for (r,_c) in colj {
            let v = self.get(r,i) + 1;
            self.set(r, i, v); 
        }
    }
    
    /// Performs the Smith normal form algorthm
    pub fn snf(&mut self) {
        self.reduce(0);
    }

    fn reduce(&mut self, i: usize) {
        // If there is some j>=i & k>=i such that m[j,k] = 1
        let search = self.entries.iter()
            .filter(|(r,c)| *r >= i && *c >= i)
            .next();

        if let Some((j,k)) = search {
            let j = *j;
            let k = *k;
            self.swap_rows(j, i);
            self.swap_columns(k, i);

            // Zero out 1's below (i,i)
            let rows: Vec<usize> = self.entries.iter()
                .filter(|(r,c)| *c == i && *r > i)
                .map(|(r,_c)| *r)
                .collect();
            for r in rows {
                self.add_row(r, i);
            }
    
            // Zero out 1's right of (i,i)
            let cols: Vec<usize> = self.entries.iter()
                .filter(|(r,c)| *r == i && *c > i)
                .map(|(_r,c)| *c)
                .collect();
            for c in cols {
                self.add_column(c, i);
            }

            // Repeat for next block
            self.reduce(i+1);
        }
    }
}

impl fmt::Display for SparseBoundaryMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let shape = self.shape();

        for i in 0..shape.0 {
            writeln!(f, "")?;
            for j in 0..shape.1 {
                write!(f, " {}", self.get(i,j))?;
            }
        }

        writeln!(f, "")
    }
}

