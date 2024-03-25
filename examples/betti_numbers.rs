extern crate scomplex;
use scomplex::*;
use std::time::Instant;

fn compute(name: &str, complex: &Complex) {
    let now = Instant::now();
    print!("{:<15}", name);
    let mut euler_characteristic = 0;
    for p in 0..=2 {
        let betti = complex.betti(p); 
        euler_characteristic += (if p % 2 == 0 { 1 } else { -1 }) * betti;
        print!(" {:3} ", betti);
    }
    let elapsed = now.elapsed();
    println!(" {:2}  {:>9?}", euler_characteristic, elapsed);
}

fn main() {
    println!("Manifold         β₀   β₁   β₂   Χ   Time");
    println!("--------------- ---- ---- ---- --- ---------");
    compute("Disk", &Complex::disk());
    compute("Sphere", &Complex::sphere());
    compute("Cylinder", &Complex::cylinder());
    compute("Mobius Strip", &Complex::mobius_strip());
    compute("Torus", &Complex::torus());
    compute("Klein Bottle", &Complex::klein_bottle());
    compute("ℝP²", &Complex::real_projective_plane());
}
