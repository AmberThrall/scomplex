extern crate scomplex;
use scomplex::*;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::time::Instant;

fn main() {
    // Load in the data
    let file = File::open("examples/triangles.txt").expect("failed to open triangles file");
    let reader = BufReader::new(file);
    let mut triangles = Vec::new();

    for line in reader.lines() {
        let indices: Vec<usize> = line.unwrap()
            .split_whitespace()
            .map(|x| x.parse::<usize>().expect("failed to parse data"))
            .collect();
        triangles.push(indices);
    }

    // Construction of complex
    let now = Instant::now();
    let mut complex = Complex::new();
    for tri in triangles {
        complex.push(Simplex::from(tri));
    }
    let elapsed = now.elapsed();
    println!("Construction: {:.2?}", elapsed);

    // Betti numbers
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
