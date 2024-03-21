use super::{
    Complex, Chain, complex::SimplexHandle,
};
use nalgebra as na;

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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Complex, splx, Chain};

    fn print_chain(chain: &Chain, complex: &Complex) {
        if chain.is_zero() { print!("0"); }
        for (i, h) in chain.0.iter().enumerate() {
            if i > 0 { print!("+") } 
            print!("{}", complex.get_simplex(*h).unwrap());
        }
    }

    #[test]
    fn boundary_matrix() {
        let mut complex = Complex::new();
        complex.push(splx![0,1]);
        complex.push(splx![0,2]);
        complex.push(splx![1,2]);
        complex.push(splx![1,3]);
        complex.push(splx![1,4]);
        complex.push(splx![2,3]);
        complex.push(splx![3,4]);
        complex.push(splx![0,1,2]);
        complex.push(splx![1,3,4]);
        
        let soln = na::DMatrix::from_row_slice(5, 7, &[
            1, 1, 0, 0, 0, 0, 0,
            1, 0, 1, 1, 1, 0, 0,
            0, 1, 1, 0, 0, 1, 0,
            0, 0, 0, 1, 0, 1, 1,
            0, 0, 0, 0, 1, 0, 1,
        ]);

        let bmatrix = BoundaryMatrix::new(1, &complex).unwrap();
        println!("Solution: {}", soln);
        println!("Got: {}", bmatrix.m);
    }
}
