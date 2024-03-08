extern crate scomplex;
use scomplex::*;
use std::time::Instant;

fn main() {
    // Create the data
    let data = vec![
        [1.0, 1.0],
        [7.0, 0.0],
        [4.0, 6.0],
        [9.0, 6.0],
        [0.0, 14.0],
        [2.0, 19.0],
        [9.0, 17.0],
    ];

    // Create the complex
    let now = Instant::now();
    let complex = AlphaComplex::delaunay(&data).build();
    let elapsed = now.elapsed();

    println!("Construction time: {:.2?}", elapsed);

    // Print the simplices
    let d = complex.dim();
    let n_simplices = complex.num_simplices();
    let n_vertices = complex.num_simplices_by_dimension(0);
    println!("Alpha complex is of dimension {} - {} simplices - {} vertices.", d, n_simplices, n_vertices);
    for k in 0..d+1 {
        for simplex in complex.iter().map(|h| complex.get(h).unwrap()).filter(|s| s.simplex.dim() == k) {
            println!("  {} -> [{}]", simplex.simplex, simplex.filtration_value);
        }
    }
}
