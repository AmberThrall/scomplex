extern crate scomplex;
use scomplex::*;
use std::time::Instant;

fn euclidean_distance_squared(a: &Point, b: &Point) -> f32 {
    let mut s: f32 = 0.0;

    for i in 0..a.len() {
        s += (a[i] - b[i]) * (a[i] - b[i]);
    }

    s
}

fn main() {
    // Create the data
    let data = vec![
        vec![-0.7269727879799237, -0.0845683083484674],
        vec![-0.30595887151355927, -0.24130148115124084],
        vec![0.8369789928055023, -0.9287987321576692],
        vec![0.8370428241258603, -0.1459067626757009],
        vec![-0.2890776075999093, -0.5621479124062236],
        vec![0.005118966838890016, 0.01637336792526245],
        vec![0.4502422267315644, -0.0022429236529895036],
        vec![0.17769137403483803, -0.7774016076589387],
        vec![0.29284937851925453, 0.49100126982128445],
        vec![-0.5577122023962757, 0.20472669982459157]
    ];
    let radius = 0.282;

    // Create the complex
    let now = Instant::now();
    let complex = RipsComplex::new(&data, 4.0 * radius * radius)
        .max_dim(2)
        .distance_fn(euclidean_distance_squared)
        .build();
    let elapsed = now.elapsed();

    println!("Construction time: {:.2?}", elapsed);

    // Print the simplices
    let d = complex.dim();
    let n_simplices = complex.num_simplices();
    let n_vertices = complex.num_simplices_by_dimension(0);
    println!("Rips complex is of dimension {} - {} simplices - {} vertices.", d, n_simplices, n_vertices);
    for k in 0..d+1 {
        for simplex in complex.iter().map(|h| complex.get(h).unwrap()).filter(|s| s.simplex.dim() == k) {
            println!("  {} -> [{}]", simplex.simplex, simplex.filtration_value);
        }
    }
}
