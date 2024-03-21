extern crate scomplex;
use scomplex::*;
use std::time::Instant;

fn main() {
    let mut complex = Complex::new();
    complex.push(splx![0,1,2]);
    complex.push(splx![1,3,4]);
    complex.push(splx![2,3]);

    let mut euler_characteristic = 0;
    for p in 0..=complex.dim() {
        let now = Instant::now();
        let betti = complex.betti(p);
        euler_characteristic += (if p % 2 == 0 { 1 } else { -1 }) * betti;
        let elapsed = now.elapsed();
        println!("beta_{} = {} (computed in {:.2?})", p, betti, elapsed);
    }

    println!("Sanity: {} == {}", euler_characteristic, complex.euler_characteristic());
}
