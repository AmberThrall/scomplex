use super::{
    Complex, Chain, complex::SimplexHandle,
};
use nalgebra as na;

/// Represents a boundary matrix of Z_2
pub struct BoundaryMatrix {
    pub col_labels: Vec<Chain>,
    pub row_labels: Vec<Chain>,
    pub m: na::DMatrix<i32>,
}

impl BoundaryMatrix {
    pub fn new(p: i32, complex: &Complex) -> Option<Self> {
        if p < 0 { return None; }

        // Get the p-simplices and (p-1)-simplices
        let p_simplices: Vec<SimplexHandle> = complex.iter()
            .filter(|h| complex[*h].simplex.dim() == p).collect();
        let pm1_simplices: Vec<SimplexHandle> = complex.iter()
            .filter(|h| complex[*h].simplex.dim() == p-1).collect();
        
        // Label rows and columns by chains.
        let col_labels: Vec<Chain> = p_simplices.iter()
            .map(|h| Chain::zero() + *h).collect();
        let row_labels: Vec<Chain> = pm1_simplices.iter()
            .map(|h| Chain::zero() + *h).collect();

        // Construct the boundary matrix.
        let m = na::DMatrix::from_fn(row_labels.len(), col_labels.len(), |i, j| {
            let sigma = complex.get_simplex(p_simplices[j]).unwrap();
            let tau = complex.get_simplex(pm1_simplices[i]).unwrap();
            if tau.is_face(&sigma) { 1 } else { 0 }
        });

        Some(BoundaryMatrix {
            col_labels, row_labels,
            m
        })
    }

    pub fn swap_columns(&mut self, i: usize, j: usize) {
        // Swap the matrix columns
        let coli = self.m.column(i).clone_owned();
        let colj = self.m.column(j).clone_owned();
        self.m.set_column(i, &colj);
        self.m.set_column(j, &coli);

        // Swap the labels
        let coli = self.col_labels[i].clone();
        let colj = self.col_labels[j].clone();
        self.col_labels[i] = colj;
        self.col_labels[j] = coli;
    }

    pub fn swap_rows(&mut self, i: usize, j: usize) {
        // Swap the matrix columns
        let rowi = self.m.row(i).clone_owned();
        let rowj = self.m.row(j).clone_owned();
        self.m.set_row(i, &rowj);
        self.m.set_row(j, &rowi);

        // Swap the labels
        let rowi = self.row_labels[i].clone();
        let rowj = self.row_labels[j].clone();
        self.row_labels[i] = rowj;
        self.row_labels[j] = rowi;
    }

    /// Performs column replacement, C_i <- C_i + C_j
    pub fn add_column(&mut self, i: usize, j: usize) {
        // Add the matrix columns
        let colj = self.m.column(j).clone_owned(); 
        let mut coli = self.m.column_mut(i); 

        for k in 0..coli.len() { 
            let v = if coli[k] == colj[k] { 0 } else { coli[k] + colj[k] };
            coli[k] = v;        
        }

        // Add the labels
        let colj = self.col_labels[j].clone();
        self.col_labels[i] += colj;
    }

    /// Performs row replacement, R_i <- R_i + R_j
    pub fn add_row(&mut self, i: usize, j: usize) {
        // Add the matrix rows
        let rowj = self.m.row(j).clone_owned(); 
        let mut rowi = self.m.row_mut(i); 

        for k in 0..rowi.len() { 
            let v = if rowi[k] == rowj[k] { 0 } else { rowi[k] + rowj[k] };
            rowi[k] = v;        
        }

        // Add the labels
        let rowj = self.row_labels[j].clone();
        self.row_labels[i] += rowj;
    }

    /// Performs the Smith normal form algorthm
    pub fn snf(&mut self) {
        self.reduce(0);
    }

    fn reduce(&mut self, i: usize) {
        let m = self.row_labels.len();
        let n = self.col_labels.len();

        // If there is some j>=i and k>=i such that m[j,k] == 1
        let mut found = false;
        'outer: for j in i..m {
            for k in i..n {
                if self.m[(j,k)] == 1 {
                    // swap rows i & j and columns i & k
                    found = true;
                    self.swap_rows(i, j);
                    self.swap_columns(i, k);
                    break 'outer;
                }
            }
        }

        if found {
            // Zero out 1's below m[i,i]
            for h in (i+1)..m {
                if self.m[(h,i)] == 1 {
                    self.add_row(h, i);
                }
            }
            
            // Zero out 1's to the right of m[i,i]
            for h in (i+1)..n {
                if self.m[(i,h)] == 1 {
                    self.add_column(h, i);
                }
            }

            // Proceed to next block
            self.reduce(i+1);
        }
    }
}
