extern crate scomplex;
use scomplex::{Complex, Simplex};
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
        complex.push_recursive(Simplex::from(tri));
    }
    let elapsed = now.elapsed();
    println!("Construction: {:.2?}", elapsed);

    // Euler Characteristic
    let now = Instant::now();
    let euler = complex.euler_characteristic();
    let elapsed = now.elapsed();
    println!("Euler Characteristic: {} ({:.2?})", euler, elapsed);

    // Orientable
    let now = Instant::now();
    let orientable = complex.orient().is_ok();
    let elapsed = now.elapsed();
    println!("Orientable: {} ({:.2?})", orientable, elapsed);
}