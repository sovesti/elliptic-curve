use std::ops::Shr;

use crypto_bigint::Uint;

pub trait MulMod<Rhs = Self> {
    type Output;

    fn mul_mod(&self, rhs: &Self, modulus: &Self) -> Self::Output;
}

impl<const LIMBS: usize> MulMod for Uint<LIMBS> {
    type Output = Self;

    fn mul_mod(&self, rhs: &Self, modulus: &Self) -> Self::Output {
        if *rhs == Uint::ZERO {
            Uint::ZERO
        } else {
            let half = rhs.shr(1);
            let is_odd = rhs.trailing_zeros();
            let doubled = self.add_mod(self, modulus);
            match is_odd {
                0 => doubled.mul_mod(&half, modulus).add_mod(self, modulus),
                _ => doubled.mul_mod(&half, modulus)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{Uint, U256};

    use crate::arithmetics::mul_mod::MulMod;

    #[test]
    fn mul_by_inv() {
        let p: U256 = Uint::from_be_hex("8000000000000000000000000000000000000000000000000000000000000431");
        let a = Uint::from_u32(2);
        assert_eq!(Uint::ONE, a.mul_mod(&a.inv_odd_mod(&p).0, &p));
    }
}