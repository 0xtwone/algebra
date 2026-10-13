use crate::models::short_weierstrass::SWCurveConfig;
use ark_ff::{Field, One, Zero};
use core::marker::PhantomData;

use crate::{
    hashing::{curve_maps::parity, map_to_curve_hasher::MapToCurve, HashToCurveError},
    models::short_weierstrass::{Affine, Projective},
};

/// Trait defining the necessary parameters for the SWU hash-to-curve method
/// for the curves of Weierstrass form of:
/// y^2 = x^3 + a*x + b where ab != 0. From [\[WB2019\]]
///
/// - [\[WB2019\]] <https://eprint.iacr.org/2019/403>
pub trait SWUConfig: SWCurveConfig {
    /// An element of the base field that is not a square root see \[WB2019, Section 4\].
    /// It is also convenient to have $g(b/ZETA * a)$ to be square. In general
    /// we use a `ZETA` with low absolute value coefficients when they are
    /// represented as integers.
    const ZETA: Self::BaseField;
}

/// Represents the SWU hash-to-curve map defined by `P`.
pub struct SWUMap<P: SWUConfig>(PhantomData<fn() -> P>);

impl<P: SWUConfig> MapToCurve<Projective<P>> for SWUMap<P> {
    /// Checks if `P` represents a valid map.
    fn check_parameters() -> Result<(), HashToCurveError> {
        // Verifying that ZETA is a non-square
        debug_assert!(
            P::ZETA.legendre().is_qnr(),
            "ZETA should be a quadratic non-residue for the SWU map"
        );

        // Verifying the prerequisite for applicability  of SWU map
        debug_assert!(!P::COEFF_A.is_zero() && !P::COEFF_B.is_zero(),
		      "Simplified SWU requires a * b != 0 in the short Weierstrass form of y^2 = x^3 + a*x + b ");

        Ok(())
    }

    /// Map an arbitrary base field element to a curve point.
    /// Based on
    /// <https://github.com/zcash/pasta_curves/blob/main/src/hashtocurve.rs>.
    fn map_to_curve(element: P::BaseField) -> Result<Affine<P>, HashToCurveError> {
        // 1. tv1 = inv0(Z^2 * u^4 + Z * u^2)
        // 2. x1 = (-B / A) * (1 + tv1)
        // 3. If tv1 == 0, set x1 = B / (Z * A)
        // 4. gx1 = x1^3 + A * x1 + B
        //
        // We use the "Avoiding inversions" optimization in [WB2019, section 4.2]
        // (not to be confused with section 4.3):
        //
        //   here       [WB2019]
        //   -------    ---------------------------------
        //   Z          ξ
        //   u          t
        //   Z * u^2    ξ * t^2 (called u, confusingly)
        //   x1         X_0(t)
        //   x2         X_1(t)
        //   gx1        g(X_0(t))
        //   gx2        g(X_1(t))
        //
        // Using the "here" names:
        //    x1 = num_x1/div      = [B*(Z^2 * u^4 + Z * u^2 + 1)] / [-A*(Z^2 * u^4 + Z * u^2]
        //   gx1 = num_gx1/div_gx1 = [num_x1^3 + A * num_x1 * div^2 + B * div^3] / div^3
        let a = P::COEFF_A;
        let b = P::COEFF_B;

        let zeta_u2 = P::ZETA * element.square();
        let ta = zeta_u2.square() + zeta_u2;
        let num_x1 = b * (ta + <P::BaseField as One>::one());
        let div = a * if ta.is_zero() { P::ZETA } else { -ta };

        let num2_x1 = num_x1.square();
        let div2 = div.square();
        let div3 = div2 * div;
        let num_gx1 = (num2_x1 + a * div2) * num_x1 + b * div3;

        // 5. x2 = Z * u^2 * x1
        let num_x2 = zeta_u2 * num_x1; // same div

        // 6. gx2 = x2^3 + A * x2 + B  [optimized out; see below]
        // 7. If is_square(gx1), set x = x1 and y = sqrt(gx1)
        // 8. Else set x = x2 and y = sqrt(gx2)
        let gx1_square;
        let gx1;

        debug_assert!(
            !div3.is_zero(),
            "we have checked that neither a or ZETA are zero. Q.E.D."
        );
        let y1: P::BaseField = {
            gx1 = num_gx1 / div3;
            if gx1.legendre().is_qr() {
                gx1_square = true;
                gx1.sqrt()
                    .expect("We have checked that gx1 is a quadratic residue. Q.E.D")
            } else {
                let zeta_gx1 = P::ZETA * gx1;
                gx1_square = false;
                zeta_gx1.sqrt().expect(
                    "ZETA * gx1 is a quadratic residue because legard is multiplicative. Q.E.D",
                )
            }
        };

        // This optimization also comes from a generalization of [WB2019, section 4.2].
        //
        // We use the specialization with h = Z = ZETA (a fixed quadratic non-residue).
        // Since gx2 = g(Z * u^2 * x1) = Z^3 * u^6 * gx1, when gx1 is not square we take
        // y1 such that y1^2 = Z * gx1 (i.e., y1 = sqrt(Z * gx1)), and then set
        // y2 = Z * u^3 * y1. This gives y2^2 = (Z * u^3)^2 * (Z * gx1) = Z^3 * u^6 * gx1 = gx2,
        // so we avoid computing gx2 explicitly.

        let y2 = zeta_u2 * element * y1;
        let num_x = if gx1_square { num_x1 } else { num_x2 };
        let y = if gx1_square { y1 } else { y2 };

        let x_affine = num_x / div;
        let y_affine = if parity(&y) == parity(&element) {
            y
        } else {
            -y
        };
        let point_on_curve = Affine::new_unchecked(x_affine, y_affine);
        debug_assert!(
            point_on_curve.is_on_curve(),
            "swu mapped to a point off the curve"
        );
        Ok(point_on_curve)
    }
}
