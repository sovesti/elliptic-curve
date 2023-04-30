use core::fmt;

use crypto_bigint::Uint;

use crate::arithmetics::{pow_mod::PowMod, mul_mod::MulMod};

use super::point::{Point, AffinePoint};

pub enum CalculusError {
    UndefinedCurve,
}

impl fmt::Display for CalculusError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            CalculusError::UndefinedCurve => write!(f, "bad raw data, curve is undefined")
        }
    }
}

fn is_defined_by_coeffs<const LIMBS: usize>(modulus: &Uint<LIMBS>, coeff_0: &Uint<LIMBS>, coeff_1: &Uint<LIMBS>) 
-> (Uint<LIMBS>, bool) {
    let d = coeff_1
        .pow_mod(&Uint::from_u32(3), modulus)
        .mul_mod(&Uint::from_u32(4), modulus)
        .add_mod(&coeff_0
            .pow_mod(&Uint::from_u32(2), modulus)
            .mul_mod(&Uint::from_u32(27), modulus), modulus);
    if d == Uint::ZERO {
        (d, false)
    } else {
        (d, true)
    }
}

fn j_invariant<const LIMBS: usize>(modulus: &Uint<LIMBS>, d: &Uint<LIMBS>, coeff_1: &Uint<LIMBS>) -> Uint<LIMBS> {
    Uint::from_u32(1728)
    .mul_mod(&coeff_1
        .pow_mod(&Uint::from_u32(3), modulus)
        .mul_mod(&Uint::from_u32(4), modulus), modulus)
    .mul_mod(&d.inv_odd_mod(modulus).0, modulus)
}

fn is_defined_by_invariant<const LIMBS: usize>(modulus: &Uint<LIMBS>, j_invariant: &Uint<LIMBS>)
-> (Uint<LIMBS>, bool) {
    if *j_invariant == Uint::ZERO || *j_invariant == Uint::from_u32(1728) {
        (Uint::ZERO, false)
    } else {
        (j_invariant
            .mul_mod(&Uint::inv_odd_mod(
                &Uint::from_u32(1728)
                .sub_mod(j_invariant, modulus), modulus).0, modulus), 
        true)    
    }
}

pub struct Curve<const LIMBS: usize> {
    modulus: Uint<LIMBS>,
    coeff_0: Uint<LIMBS>,
    coeff_1: Uint<LIMBS>,
    j_invariant: Uint<LIMBS>
}

impl<const LIMBS: usize> Curve<LIMBS> {
    pub fn from_coeffs(modulus: Uint<LIMBS>, coeff_0: Uint<LIMBS>, coeff_1: Uint<LIMBS>) ->
    Result<Curve<LIMBS>, CalculusError> {
        let (d, is_defined) = is_defined_by_coeffs(&modulus, &coeff_0, &coeff_1);
        if is_defined {
            Ok(Curve {
                modulus: modulus,
                coeff_0: coeff_0,
                coeff_1: coeff_1,
                j_invariant: j_invariant(&modulus, &d, &coeff_1)
            })
        } else {
            Err(CalculusError::UndefinedCurve)
        }
    }

    pub fn from_invariant(modulus: Uint<LIMBS>, j_invariant: Uint<LIMBS>) ->
    Result<Curve<LIMBS>, CalculusError> {
        let (k, is_defined) = is_defined_by_invariant(&modulus, &j_invariant);
        if is_defined {
            Ok(Curve { 
                modulus: modulus, 
                coeff_0: k.mul_mod(&Uint::from_u32(2), &modulus), 
                coeff_1: k.mul_mod(&Uint::from_u32(3), &modulus), 
                j_invariant: j_invariant 
            })
        } else {
            Err(CalculusError::UndefinedCurve)
        }
    }

    pub fn contains_point(&self, point: Point<LIMBS>) -> bool {
        match point {
            Point::Affine(point) => 
                point.y.mul_mod(&point.y, &self.modulus) == 
                point.x.pow_mod(&Uint::from_u8(3), &self.modulus)
                .add_mod(&point.x.mul_mod(&self.coeff_1, &self.modulus), &self.modulus)
                .add_mod(&self.coeff_0, &self.modulus),
            Point::Projective => true,
        }
    }

    pub fn add_points(&self, lhs: Point<LIMBS>, rhs: Point<LIMBS>) -> Point<LIMBS> {
        match lhs {
            Point::Affine(lhs) => {
                match rhs {
                    Point::Affine(rhs) => self.add_affine_points(lhs, rhs),
                    Point::Projective => Point::Affine(lhs),
                }
            },
            Point::Projective => rhs,
        }
    }

    fn sum_from_lambda(&self, lhs: AffinePoint<LIMBS>, rhs: AffinePoint<LIMBS>, lambda: Uint<LIMBS>)
    -> Point<LIMBS> {
        let result_x = lambda
            .mul_mod(&lambda, &self.modulus)
            .sub_mod(&lhs.x, &self.modulus)
            .sub_mod(&rhs.x, &self.modulus);
        let result_y = lambda
        .mul_mod(&lhs.x
            .sub_mod(&result_x, &self.modulus), &self.modulus)
        .sub_mod(&lhs.y, &self.modulus);
        Point::Affine(AffinePoint::new(result_x, result_y))
    }

    fn add_affine_points(&self, lhs: AffinePoint<LIMBS>, rhs: AffinePoint<LIMBS>) 
    -> Point<LIMBS> {
        if lhs.x == rhs.x {
            if lhs.y == Uint::ZERO.sub_mod(&rhs.y, &self.modulus) {
                Point::Projective
            } else {
                let lambda = lhs.x
                .mul_mod(&lhs.x, &self.modulus)
                .mul_mod(&Uint::from_u32(3), &self.modulus)
                .add_mod(&self.coeff_1, &self.modulus)
                .mul_mod(&lhs.y
                    .mul_mod(&Uint::from_u32(2), &self.modulus)
                    .inv_odd_mod(&self.modulus).0, &self.modulus);
                self.sum_from_lambda(lhs, rhs, lambda)
            }
        } else {
            let lambda = rhs.y
            .sub_mod(&lhs.y, &self.modulus)
            .mul_mod(&rhs.x
                .sub_mod(&lhs.x, &self.modulus)
                .inv_odd_mod(&self.modulus).0, &self.modulus);
            self.sum_from_lambda(lhs, rhs, lambda)
        }
    }
}

