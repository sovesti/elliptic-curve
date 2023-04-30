use crypto_bigint::Uint;
use std::ops::Shr;
use crate::arithmetics::mul_mod::MulMod;

pub trait PowMod<Rhs = Self> {
    type Output;

    fn pow_mod(&self, rhs: &Rhs, modulus: &Self) -> Self::Output;
}

impl<const LIMBS: usize> PowMod for Uint<LIMBS> {
    type Output = Self;

    fn pow_mod(&self, rhs: &Self, modulus: &Self) -> Self::Output {
        if *rhs == Uint::ZERO {
            Uint::ONE
        } else {
            let half = rhs.shr(1);
            let is_odd = rhs.trailing_zeros();
            let squared = self.mul_mod(self, modulus);
            match is_odd {
                0 => squared.pow_mod(&half, modulus).mul_mod(self, modulus),
                _ => squared.pow_mod(&half, modulus)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crypto_bigint::{Uint, U256};

    use crate::arithmetics::pow_mod::PowMod;

    #[test]
    fn fermat_little_theorem() {
        let p: U256 = Uint::from_be_hex("8000000000000000000000000000000000000000000000000000000000000431");
        let a = Uint::from_u32(2);
        assert_eq!(a.pow_mod(&p.sub_mod(&Uint::ONE, &p), &p), Uint::ONE);
    }
}