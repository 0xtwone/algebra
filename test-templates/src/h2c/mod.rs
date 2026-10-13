pub mod json;
use hex;
pub use hex::decode;
use serde_json;
pub use serde_json::from_reader;
use sha2;
pub use sha2::Sha256;

/// Checks a hash-to-curve suite against the known-answer test vectors of
/// [RFC 9380](https://datatracker.ietf.org/doc/html/rfc9380), read from
/// `{test_path}/{test_name}_XMD-SHA-256_SSWU_RO_.json`.
///
/// Test vectors only exist for a handful of curves; for randomized tests that
/// work with any curve, see [`test_h2c!`](crate::test_h2c).
#[macro_export]
macro_rules! test_h2c_suite {
    ($mod_name: ident; $test_path: literal; $test_name: literal; $group: ty; $field: ty; $base_prime_field: ty; $m: literal) => {
        mod $mod_name {
            use ark_ff::PrimeField;

            extern crate std;
            use ark_ec::{
                hashing::{
                    curve_maps::wb::WBMap, map_to_curve_hasher::MapToCurveBasedHasher, HashToCurve,
                },
                short_weierstrass::{Affine, Projective},
            };
            use ark_ff::{
                field_hashers::{DefaultFieldHasher, HashToField},
                fields::Field,
                One, UniformRand, Zero,
            };
            use ark_std::{format, string::String, vec::*};
            use std::{
                fs::{read_dir, File},
                io::BufReader,
            };
            use $crate::{decode, Sha256};

            use $crate::json::SuiteVector;
            #[test]
            fn test_h2c() {
                let filename = format!("{}/{}_XMD-SHA-256_SSWU_RO_.json", $test_path, $test_name);

                let file = File::open(filename).unwrap();
                let data: SuiteVector = $crate::from_reader(BufReader::new(file)).unwrap();

                assert_eq!(data.hash, "sha256");
                let dst = data.dst.as_bytes();
                let hasher;
                let g1_mapper = MapToCurveBasedHasher::<
                    Projective<$group>,
                    DefaultFieldHasher<Sha256, 128>,
                    WBMap<$group>,
                >::new(dst)
                .unwrap();
                hasher = <DefaultFieldHasher<Sha256, 128> as HashToField<$field>>::new(dst);

                for v in data.vectors.iter() {
                    // first, hash-to-field tests
                    let got: [$base_prime_field; { 2 * $m }] =
                        hasher.hash_to_field(&v.msg.as_bytes());
                    let want: Vec<$base_prime_field> =
                        v.u.iter().map(read_fq_vec).flatten().collect();
                    assert_eq!(got[..], *want);

                    // then, test curve points
                    let x = read_fq_vec(&v.p.x);
                    let y = read_fq_vec(&v.p.y);
                    let got = g1_mapper.hash(&v.msg.as_bytes()).unwrap();
                    let want = Affine::<$group>::new_unchecked(
                        <$field>::from_base_prime_field_elems(x).unwrap(),
                        <$field>::from_base_prime_field_elems(y).unwrap(),
                    );
                    assert!(got.is_on_curve());
                    assert!(want.is_on_curve());
                    assert_eq!(got, want);
                }
            }
            pub fn read_fq_vec(input: &String) -> Vec<$base_prime_field> {
                input
                    .split(",")
                    .map(|f| {
                        <$base_prime_field>::from_be_bytes_mod_order(
                            &decode(f.trim_start_matches("0x")).unwrap(),
                        )
                    })
                    .collect()
            }
        }
    };
}

/// Randomized tests for a [`MapToCurve`] implementation and the
/// [`MapToCurveBasedHasher`] built on top of it. Unlike
/// [`test_h2c_suite!`](crate::test_h2c_suite), these need no test vectors, so
/// they can be run for any curve that supports hashing.
///
/// Usage: `test_h2c!(mod_name; Config; MapType; model)`, where `Config` is
/// the configuration of the curve that `MapType` maps to and `model` is the
/// curve model module in `ark_ec::models`, i.e. `short_weierstrass` or
/// `twisted_edwards`.
///
/// [`MapToCurve`]: ark_ec::hashing::map_to_curve_hasher::MapToCurve
/// [`MapToCurveBasedHasher`]: ark_ec::hashing::map_to_curve_hasher::MapToCurveBasedHasher
#[macro_export]
macro_rules! test_h2c {
    ($mod_name: ident; $config: ty; $map: ty; $model: ident) => {
        mod $mod_name {
            extern crate std;
            use ark_ec::{
                hashing::{
                    map_to_curve_hasher::{MapToCurve, MapToCurveBasedHasher},
                    HashToCurve,
                },
                $model::{Affine, Projective},
                CurveConfig,
            };
            use ark_ff::{field_hashers::DefaultFieldHasher, UniformRand};
            use ark_std::{rand::Rng, test_rng, vec::Vec};
            use std::collections::HashSet;
            use $crate::Sha256;

            type BaseField = <$config as CurveConfig>::BaseField;
            type Hasher =
                MapToCurveBasedHasher<Projective<$config>, DefaultFieldHasher<Sha256, 128>, $map>;

            const ITERATIONS: usize = 100;
            const DST: &[u8] = b"ARKWORKS-H2C-TEST-V01";

            fn map_to_curve(u: BaseField) -> Affine<$config> {
                <$map as MapToCurve<Projective<$config>>>::map_to_curve(u).unwrap()
            }

            fn random_message<R: Rng>(rng: &mut R) -> Vec<u8> {
                let len = rng.gen_range(0..256);
                (0..len).map(|_| rng.gen()).collect()
            }

            #[test]
            fn check_parameters() {
                <$map as MapToCurve<Projective<$config>>>::check_parameters().unwrap();
            }

            #[test]
            fn map_to_curve_is_on_curve() {
                let mut rng = test_rng();
                for _ in 0..ITERATIONS {
                    let u = BaseField::rand(&mut rng);
                    let point = map_to_curve(u);
                    assert!(point.is_on_curve());
                    assert_eq!(point, map_to_curve(u));
                }
            }

            #[test]
            fn map_to_curve_is_not_constant() {
                let mut rng = test_rng();
                let images: HashSet<_> = (0..ITERATIONS)
                    .map(|_| map_to_curve(BaseField::rand(&mut rng)))
                    .collect();
                assert!(images.len() > ITERATIONS / 2);
            }

            #[test]
            fn hash_to_curve_is_in_correct_subgroup() {
                let mut rng = test_rng();
                let hasher = Hasher::new(DST).unwrap();
                for _ in 0..ITERATIONS {
                    let msg = random_message(&mut rng);
                    let point = hasher.hash(&msg).unwrap();
                    assert!(point.is_on_curve());
                    assert!(point.is_in_correct_subgroup_assuming_on_curve());
                    assert_eq!(point, hasher.hash(&msg).unwrap());
                }
            }

            #[test]
            fn hash_to_curve_depends_on_message_and_dst() {
                let mut rng = test_rng();
                let hasher = Hasher::new(DST).unwrap();
                let other_hasher = Hasher::new(b"ARKWORKS-H2C-TEST-V02").unwrap();
                let mut images = HashSet::new();
                for i in 0..ITERATIONS {
                    // Prefix with the counter so that no two messages are equal.
                    let mut msg = (i as u64).to_le_bytes().to_vec();
                    msg.extend(random_message(&mut rng));
                    let point = hasher.hash(&msg).unwrap();
                    assert_ne!(point, other_hasher.hash(&msg).unwrap());
                    assert!(images.insert(point));
                }
            }
        }
    };
}
