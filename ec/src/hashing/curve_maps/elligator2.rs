use crate::models::twisted_edwards::{MontCurveConfig, TECurveConfig};
use ark_ff::{Field, One, Zero};
use core::marker::PhantomData;

use crate::{
    hashing::{curve_maps::parity, map_to_curve_hasher::MapToCurve, HashToCurveError},
    models::twisted_edwards::{Affine, Projective},
};

/// Trait defining the necessary parameters for the Elligator2 hash-to-curve method
/// for twisted edwards curves form of:
/// `b * y² = x³ + a * x² + x`
/// from [\[BHKL13\]], according to [\[HSSWW23\]]
///
/// - [\[BHKL13\]] <http://dx.doi.org/10.1145/2508859.2516734>
/// - [\[HSSWW23\]] <https://datatracker.ietf.org/doc/html/rfc9380>
pub trait Elligator2Config: TECurveConfig + MontCurveConfig {
    /// An element of the base field that is not a square root see \[BHKL13, Section 5\].
    /// When `BaseField` is a prime field, [\[HSSWW23\]] mandates that `Z` is the
    /// non-square with lowest absolute value in the `BaseField` when its elements
    /// are represented as [-(q-1)/2, (q-1)/2]
    const Z: Self::BaseField;

    /// This must be equal to 1/(MontCurveConfig::COEFF_B)^2;
    const ONE_OVER_COEFF_B_SQUARE: Self::BaseField;

    /// This must be equal to MontCurveConfig::COEFF_A/MontCurveConfig::COEFF_B;
    const COEFF_A_OVER_COEFF_B: Self::BaseField;
}

/// Represents the Elligator2 hash-to-curve map defined by `P`.
pub struct Elligator2Map<P: TECurveConfig>(PhantomData<fn() -> P>);

impl<P: Elligator2Config> MapToCurve<Projective<P>> for Elligator2Map<P> {
    /// Checks if `P` represents a valid Elligator2 map. Panics otherwise.
    fn check_parameters() -> Result<(), HashToCurveError> {
        // We assume that the Montgomery curve is correct and  as such we do
        // not verify the prerequisite for applicability of Elligator2 map to the TECurveConfing.

        // Verifying that Z is a non-square
        debug_assert!(
            !P::Z.legendre().is_qr(),
            "Z should be a quadratic non-residue for the Elligator2 map"
        );

        debug_assert_eq!(
            P::ONE_OVER_COEFF_B_SQUARE,
            <P as MontCurveConfig>::COEFF_B
                .square()
                .inverse()
                .expect("B coefficient cannot be zero in Montgomery form"),
            "ONE_OVER_COEFF_B_SQUARE is not equal to 1/COEFF_B^2 in Montgomery form"
        );

        debug_assert_eq!(
            P::COEFF_A_OVER_COEFF_B,
            <P as MontCurveConfig>::COEFF_A / <P as MontCurveConfig>::COEFF_B,
            "COEFF_A_OVER_COEFF_B is not equal to COEFF_A/COEFF_B in Montgomery form"
        );
        Ok(())
    }

    /// Map an arbitrary base field element `element` to a curve point.
    fn map_to_curve(element: P::BaseField) -> Result<Affine<P>, HashToCurveError> {
        // 1. x1 = -(J / K) * inv0(1 + Z * u^2)
        // 2. If x1 == 0, set x1 = -(J / K)
        // 3. gx1 = x1^3 + (J / K) * x1^2 + x1 / K^2
        // 4. x2 = -x1 - (J / K)
        // 5. gx2 = x2^3 + (J / K) * x2^2 + x2 / K^2
        // 6. If is_square(gx1), set x = x1, y = sqrt(gx1) with sgn0(y) == 1.
        // 7. Else set x = x2, y = sqrt(gx2) with sgn0(y) == 0.
        // 8. s = x * K
        // 9. t = y * K
        // 10. return (s, t)

        // ark a is irtf J
        // ark b is irtf k
        let k = <P as MontCurveConfig>::COEFF_B;
        let j_on_k = P::COEFF_A_OVER_COEFF_B;
        let ksq_inv = P::ONE_OVER_COEFF_B_SQUARE;

        let den_1 = <P::BaseField as One>::one() + P::Z * element.square();

        let x1 = -j_on_k
            / (if den_1.is_zero() {
                <P::BaseField as One>::one()
            } else {
                den_1
            });
        let x1sq = x1.square();
        let x1cb = x1sq * x1;
        let gx1 = x1cb + j_on_k * x1sq + x1 * ksq_inv;

        let x2 = -x1 - j_on_k;
        let x2sq = x2.square();
        let x2cb = x2sq * x2;
        let gx2 = x2cb + j_on_k * x2sq + x2 * ksq_inv;

        let (x, mut y, sgn0) = if gx1.legendre().is_qr() {
            (
                x1,
                gx1.sqrt()
                    .expect("We have checked that gx1 is a quadratic residue. Q.E.D"),
                true,
            )
        } else {
            (
                x2,
                gx2.sqrt()
                    .expect("gx2 is a quadratic residue because gx1 is not. Q.E.D"),
                false,
            )
        };

        if parity(&y) != sgn0 {
            y = -y;
        }

        let s = x * k;
        let t = y * k;

        // `(s, t)` is an affine point on the Montgomery curve.
        // Ideally, the TECurve would come with a mapping to its Montgomery curve,
        // so we could just call that mapping here.
        // This is currently not supported in arkworks, so
        // we just implement the rational map here from [\[HSSWW23\]] Appendix D

        let tv1 = s + <P::BaseField as One>::one();
        let tv2 = tv1 * t;
        let (v, w) = if tv2.is_zero() {
            (<P::BaseField as Zero>::zero(), <P::BaseField as One>::one())
        } else {
            let tv2_inv = tv2
                .inverse()
                .expect("None zero element has inverse. Q.E.D.");
            let v = tv2_inv * tv1 * s;
            let w = tv2_inv * t * (s - <P::BaseField as One>::one());
            (v, w)
        };

        let point_on_curve = Affine::new_unchecked(v, w);
        debug_assert!(
            point_on_curve.is_on_curve(),
            "Elligator2 mapped to a point off the curve"
        );
        Ok(point_on_curve)
    }
}
