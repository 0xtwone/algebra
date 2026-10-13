use core::marker::PhantomData;

use crate::{models::short_weierstrass::SWCurveConfig, CurveConfig};
use ark_ff::{batch_inversion, Field};

use crate::{
    hashing::{map_to_curve_hasher::MapToCurve, HashToCurveError},
    models::short_weierstrass::{Affine, Projective},
    AffineRepr,
};

use super::swu::{SWUConfig, SWUMap};
type BaseField<MP> = <MP as CurveConfig>::BaseField;

/// [`IsogenyMap`] defines an isogeny between curves of
/// form `Phi(x, y) := (a(x), b(x)*y)`.
///
/// The `x` coordinate of the codomain point only depends on the
/// `x`-coordinate of the domain point, and the
/// `y`-coordinate of the codomain point is a multiple of the `y`-coordinate of the domain point.
/// The multiplier depends on the `x`-coordinate of the domain point.
/// All isogeny maps of curves of short Weierstrass form can be written in this way. See
/// [\[Ga18]\]. Theorem 9.7.5 for details.
///
/// We assume that `Domain` and `Codomain` have the same `BaseField` but we use both
/// `BaseField<Domain>` and `BaseField<Codomain>` in our fields' definitions to avoid
/// using `PhantomData`
///
/// - [\[Ga18]\] Galbraith, S. D. (2018). Mathematics of public key cryptography.
pub struct IsogenyMap<
    'a,
    Domain: SWCurveConfig,
    Codomain: SWCurveConfig<BaseField = BaseField<Domain>>,
> {
    pub x_map_numerator: &'a [BaseField<Domain>],
    pub x_map_denominator: &'a [BaseField<Codomain>],

    pub y_map_numerator: &'a [BaseField<Domain>],
    pub y_map_denominator: &'a [BaseField<Codomain>],
}

impl<Domain, Codomain> IsogenyMap<'_, Domain, Codomain>
where
    Domain: SWCurveConfig,
    Codomain: SWCurveConfig<BaseField = BaseField<Domain>>,
{
    fn apply(&self, domain_point: Affine<Domain>) -> Result<Affine<Codomain>, HashToCurveError> {
        match domain_point.xy() {
            Some((x, y)) => {
                let mut v: [BaseField<Domain>; 2] = [
                    horner(self.x_map_denominator, &x),
                    horner(self.y_map_denominator, &x),
                ];
                batch_inversion(&mut v);
                let img_x = horner(self.x_map_numerator, &x) * v[0];
                let img_y = (horner(self.y_map_numerator, &x) * y) * v[1];
                Ok(Affine::new_unchecked(img_x, img_y))
            },
            None => Ok(Affine::identity()),
        }
    }
}

/// Evaluates the polynomial with the given coefficients (lowest degree first)
/// at `x` by Horner's rule, without allocating.
fn horner<F: Field>(coeffs: &[F], x: &F) -> F {
    let mut iter = coeffs.iter().rev();
    let mut acc = match iter.next() {
        Some(c) => *c,
        None => return F::zero(),
    };
    for c in iter {
        acc *= x;
        acc += c;
    }
    acc
}

/// Trait defining the necessary parameters for the WB hash-to-curve method.
///
/// This method is used for curves in Weierstrass form defined by:
///
/// `y^2 = x^3 + a*x + b` where `b != 0`, but `a` can be zero,
/// as seen in curves like BLS-381.
///
/// For more information, refer to \[WB2019\].
///
/// - [\[WB2019\]] <http://dx.doi.org/10.46586/tches.v2019.i4.154-179>
pub trait WBConfig: SWCurveConfig + Sized {
    // The isogenous curve should be defined over the same base field but it can have
    // different scalar field type IsogenousCurveScalarField :
    type IsogenousCurve: SWUConfig<BaseField = BaseField<Self>>;

    const ISOGENY_MAP: IsogenyMap<'static, Self::IsogenousCurve, Self>;
}

pub struct WBMap<P: WBConfig> {
    swu_field_curve_hasher: PhantomData<SWUMap<P::IsogenousCurve>>,
    curve_params: PhantomData<fn() -> P>,
}

impl<P: WBConfig> MapToCurve<Projective<P>> for WBMap<P> {
    /// Checks if `P` represents a valid map.
    fn check_parameters() -> Result<(), HashToCurveError> {
        match P::ISOGENY_MAP.apply(P::IsogenousCurve::GENERATOR) {
            Ok(point_on_curve) => {
                debug_assert!(point_on_curve.is_on_curve(),
			      "the isogeny maps the generator of its domain: {} into {} which does not belong to its codomain.",P::IsogenousCurve::GENERATOR, point_on_curve);
            },
            Err(e) => return Err(e),
        }

        SWUMap::<P::IsogenousCurve>::check_parameters().unwrap(); // Or ?
        Ok(())
    }

    /// Map random field point to a random curve point
    /// inspired from
    /// <https://github.com/zcash/pasta_curves/blob/main/src/hashtocurve.rs>
    fn map_to_curve(
        element: <Affine<P> as AffineRepr>::BaseField,
    ) -> Result<Affine<P>, HashToCurveError> {
        // first we need to map the field point to the isogenous curve
        let point_on_isogenious_curve = SWUMap::map_to_curve(element).unwrap();
        P::ISOGENY_MAP.apply(point_on_isogenious_curve)
    }
}
