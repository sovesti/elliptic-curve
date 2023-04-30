use crypto_bigint::Uint;

pub enum Point<const LIMBS: usize> {
    Affine(AffinePoint<LIMBS>),
    Projective
}

pub struct AffinePoint<const LIMBS: usize> {
    pub x: Uint<LIMBS>,
    pub y: Uint<LIMBS>,
}

impl<const LIMBS: usize> AffinePoint<LIMBS> {
    pub fn new(x: Uint<LIMBS>, y: Uint<LIMBS>) -> Self {
        AffinePoint { x: x, y: y }        
    }
}