use crate::*;
use ark_algebra_test_templates::*;

test_group!(te; EdwardsProjective; te);
test_h2c!(elligator2; crate::BandersnatchConfig; ark_ec::hashing::curve_maps::elligator2::Elligator2Map<crate::BandersnatchConfig>; twisted_edwards);
