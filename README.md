# scomplex

Simplicial complex construction library written in Rust.

## Features
- Construction of Vietoris-Rips filtration
- Construction of Alpha filtration for 2D point clouds
- Computation of Betti numbers
- Algorithm for orienting orientable complexes

## Example

```rust
let complex = RipsComplex::new(&data, 4.0 * radius * radius)
    .max_dim(2)
    .distance_fn(euclidean_distance_squared)
    .build();

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
```
