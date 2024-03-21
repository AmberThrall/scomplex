extern crate scomplex;
use scomplex::*;

fn print_chain(chain: &Chain, complex: &Complex) {
    if chain.is_zero() { print!("0"); }
    for (i, h) in chain.0.iter().enumerate() {
        if i > 0 { print!("+") } 
        print!("{}", complex.get_simplex(*h).unwrap());
    }
}

fn print_bmatrix(bmatrix: &BoundaryMatrix, complex: &Complex) {
    print!("Col-labels:");
    for c in bmatrix.col_labels.iter() {
        print!(" ");
        print_chain(c, &complex);
    }
    println!("");

    print!("Row-labels:");
    for c in bmatrix.row_labels.iter() {
        print!(" ");
        print_chain(c, &complex);
    }
    println!("");
    println!("{}", bmatrix.m);
}

fn main() {
    let mut complex = Complex::new();
    complex.push(splx![0,1,2]);
    complex.push(splx![1,3,4]);
    complex.push(splx![2,3]);

    for p in 0..=complex.dim() {
        println!("SNF([∂{}]):", p);
        let mut bmatrix = BoundaryMatrix::new(p, &complex).unwrap();
        bmatrix.snf();
        print_bmatrix(&bmatrix, &complex);
    }
}
