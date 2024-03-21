use crate::Complex;
use crate::complex::SimplexHandle;
use std::collections::HashSet;
use std::fmt;
use std::ops;

/// Represents a chain with underlying group Z_2.
#[derive(PartialEq, Clone)]
pub struct Chain(HashSet<SimplexHandle>);

impl Chain {
    /// Constructs a zero-chain
    pub fn zero() -> Chain {
        Chain(HashSet::new())
    }

    /// Constructs a chain given a set of simplices
    pub fn new(simplices: HashSet<SimplexHandle>) -> Chain {
        Chain(simplices)
    }

    /// Checks if the chain is the zero chain.
    pub fn is_zero(&self) -> bool {
        self.0.len() == 0
    }

    /// Checks if the chain is a cycle, i.e., if the boundary is zero.
    pub fn is_cycle(&self, complex: &Complex) -> bool {
        self.boundary(complex).is_zero()
    }

    /// Finds and returns the boundary of the chain
    ///
    /// If the chain has simplices not in the complex, they are ignored.
    pub fn boundary(&self, complex: &Complex) -> Chain {
        let mut bdry = Chain::zero();
        for h in self.0.iter() {
            if let Some(sigma) = complex.get_simplex(*h) {
                if sigma.dim() == 0 { continue; } 

                for face in complex.get_parents(*h) {
                    bdry += face; 
                }
            }
        }

        bdry
    }
}

impl fmt::Debug for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Chain(")?;
        for (i, h) in self.0.iter().enumerate() {
            if i > 0 { write!(f, "+",)?; }
            write!(f, "{:?}", h)?;
        }
        write!(f, ")")
    }
}

impl ops::Add<Chain> for Chain {
    type Output = Chain;

    fn add(self, rhs: Chain) -> Chain {
        let mut clone = self.clone();
        for h in rhs.0.iter() {
            if !clone.0.insert(*h) { // chain already contains simplex (1+1 = 0)
                clone.0.remove(h); 
            }
        }
        clone
    }
}

impl ops::AddAssign<Chain> for Chain {
    fn add_assign(&mut self, rhs: Chain) {
        for h in rhs.0.iter() {
            if !self.0.insert(*h) { // chain already contains simplex (1+1 = 0)
                self.0.remove(h); 
            }
        }
    }
}

impl ops::Add<SimplexHandle> for Chain {
    type Output = Chain;

    fn add(self, rhs: SimplexHandle) -> Chain {
        let mut set = HashSet::new();
        set.insert(rhs);
        let rhs_chain = Chain::new(set);
        
        self.add(rhs_chain)
    }
}

impl ops::AddAssign<SimplexHandle> for Chain {
    fn add_assign(&mut self, rhs: SimplexHandle) {
        let mut set = HashSet::new();
        set.insert(rhs);
        let rhs_chain = Chain::new(set);
        
        self.add_assign(rhs_chain)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Complex;
    use crate::splx;

    fn print_chain(chain: &Chain, complex: &Complex) {
        if chain.is_zero() { print!("0"); }
        for (i, h) in chain.0.iter().enumerate() {
            if i > 0 { print!("+") } 
            print!("{}", complex.get_simplex(*h).unwrap());
        }
    }

    #[test]
    fn chain() {
        let mut complex = Complex::new();
        let mut chain = Chain::zero();

        chain += complex.push(splx![0,1]);
        chain += complex.push(splx![1,2]);
        chain += complex.push(splx![0,2]);

        print!("Chain: ");
        print_chain(&chain, &complex);
        println!("");

        assert!(chain.is_cycle(&complex));
    }
}

