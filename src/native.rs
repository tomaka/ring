// Copyright 2018 Pierre Krieger.
//
// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHORS DISCLAIM ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
// SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
// OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
// CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

//! This module contains native implementations of functions otherwise implemented in
//! platform-specific assembly code.

#![cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "arm", target_arch = "aarch64")))]

use std::mem;
use std::slice;
use c;
use limb::Limb;
use rsa::bigint::N0;
use bigint::uint::{U256, U512};

#[no_mangle]
pub unsafe extern "C" fn GFp_bn_mul_mont(r: *mut Limb, a: *const Limb,
                                         b: *const Limb, n: *const Limb,
                                         n0: &N0, num_limbs: c::size_t)
{
    let num_values = num_limbs * mem::size_of::<Limb>() / mem::size_of::<U256>();

    let a = slice::from_raw_parts(a as *const U256, num_values);
    let b = slice::from_raw_parts(b as *const U256, num_values);
    let n = slice::from_raw_parts(n as *const U256, num_values);
    let n0 = ((n0[0] as u64) << 32) | (n0[1] as u64);
    let n0 = U256::from(n0);
    let r = slice::from_raw_parts_mut(r as *mut U256, num_values);

    // Note that `r` can alias `a` or `b`.

    for limb in 0 .. num_values {
        let m = ((a[limb] % b[limb]).full_mul(n0)) % U512::from(b[limb]);
        let t = (U512::from(a[limb]) + m * U512::from(n[limb])) / U512::from(b[limb]);
        let out = t % U512::from(n[limb]);
        r[limb] = U256::from(out);
    }
}

const N: U256 = U256([0xffffffffffffffff, 0x00000000ffffffff, 0x0000000000000000, 0xffffffff00000001]);
const A_MONT: U256 = U256([0xfffffffffffffffc, 0x00000003ffffffff, 0x0000000000000000, 0xfffffffc00000004]);

#[no_mangle]
pub unsafe extern "C" fn GFp_nistz256_neg(r: *mut Limb, a: *const Limb) {
    let a = a as *const U256;
    let r = r as *mut U256;
    *r = neg(*a);
}

fn neg(a: U256) -> U256 {
    N.overflowing_sub(a).0 % N
}

#[no_mangle]
pub unsafe extern "C" fn GFp_nistz256_add(r: *mut Limb/*[COMMON_OPS.num_limbs]*/,
                                      a: *const Limb/*[COMMON_OPS.num_limbs]*/,
                                      b: *const Limb/*[COMMON_OPS.num_limbs]*/)
{
    let a = a as *const U256;
    let b = b as *const U256;
    let r = r as *mut U256;
    *r = add(*a, *b);
}

fn add(a: U256, b: U256) -> U256 {
	U256::from((U512::from(a) + U512::from(b)) % U512::from(N))
}

#[no_mangle]
pub unsafe extern "C" fn GFp_nistz256_mul_mont(r: *mut Limb/*[COMMON_OPS.num_limbs]*/,
                                           a: *const Limb/*[COMMON_OPS.num_limbs]*/,
                                           b: *const Limb/*[COMMON_OPS.num_limbs]*/)
{
    let a = a as *const U256;
    let b = b as *const U256;
    let r = r as *mut U256;
    *r = mul(*a, *b);
}

fn mul(a: U256, b: U256) -> U256 {
    // TODO: don't parse decimal every time
	let r_inv = U256::from_dec_str("115792089183396302114378112356516095823261736990586219612555396166510339686400").unwrap();

	let mul = (a.full_mul(b) % U512::from(N)) * U512::from(r_inv) % U512::from(N);
    mul.into()
}

#[no_mangle]
pub unsafe extern "C" fn GFp_nistz256_sqr_mont(r: *mut Limb/*[COMMON_OPS.num_limbs]*/,
                                           a: *const Limb/*[COMMON_OPS.num_limbs]*/)
{
    GFp_nistz256_mul_mont(r, a, a);
}
fn squared(a: U256) -> U256 {
	mul(a, a)
}

fn cubed(a: U256) -> U256 {
	mul(a, squared(a))
}

fn sub(a: U256, b: U256) -> U256 {
	add(a, neg(b))
}

fn doubled(a: U256) -> U256 {
	add(a, a)
}

fn quadrupled(a: U256) -> U256 {
	add(a, tripled(a))
}

fn tripled(a: U256) -> U256 {
	add(a, doubled(a))
}

fn p_double(x: U256, y: U256, z: U256) -> (U256, U256, U256) {
	if y == U256::zero() {
		return (U256::zero(), U256::zero(), U256::zero())
	}

	// S = 4*X*Y^2
	let s = mul(quadrupled(x), squared(y));

	// M = 3*X^2 + a*Z^4
	let m = add(tripled(squared(x)), mul(A_MONT, mul(cubed(z), z)));

	// X' = M^2 - 2*S
	let r_x = sub(squared(m), doubled(s));
	// Y' = M*(S - X') - 8*Y^4
	let r_y = sub(
		mul(m, sub(s, r_x)),
		doubled(quadrupled(mul(cubed(y), y))),
	);
	// Z' = 2*Y*Z
	let r_z = doubled(mul(y, z));

	(r_x, r_y, r_z)
}

#[no_mangle]
pub unsafe extern "C" fn GFp_nistz256_point_add(r: *mut Limb/*[3][COMMON_OPS.num_limbs]*/,
                                                a: *const Limb/*[3][COMMON_OPS.num_limbs]*/,
                                                b: *const Limb/*[3][COMMON_OPS.num_limbs]*/)
{
    let a = a as *const [U256; 3];
    let b = b as *const [U256; 3];
    let r = r as *mut [U256; 3];

    let a = *a;
    let b = *b;

    let out = p_add(a[0], a[1], a[2], b[0], b[1], b[2]);
    *r = [out.0, out.1, out.2];
}

fn p_add(a_x: U256, a_y: U256, a_z: U256, b_x: U256, b_y: U256, b_z: U256) -> (U256, U256, U256) {
    if b_z.is_zero() || a_z.is_zero() {
         if b_z.is_zero() && a_z.is_zero() {
             return (0.into(), 0.into(), 0.into())
         } else if b_z.is_zero() {
             return (a_x, a_y, a_z)
         } else if a_z.is_zero() {
             return (b_x, b_y, b_z)
         }
     }

	// U1 = X1*Z2^2
	let u1 = mul(a_x, squared(b_z));
	// U2 = X2*Z1^2
	let u2 = mul(b_x, squared(a_z));

	// S1 = Y1*Z2^3
	let s1 = mul(a_y, cubed(b_z));
	// S2 = Y2*Z1^3
	let s2 = mul(b_y, cubed(a_z));

	if u1 == u2 {
		if s1 != s2 {
			return (U256::zero(), U256::zero(), U256::zero())
		} else {
			return p_double(a_x, a_y, a_z);
		}
	}

	// H = U2 - U1
	let h = sub(u2, u1);
	// R = S2 - S1
	let r = sub(s2, s1);

	// X3 = R^2 - H^3 - 2*U1*H^2
	let dd = mul(u1, squared(h));
	let r_x = sub(sub(squared(r), cubed(h)), doubled(dd));

    // Y3 = R*(U1*H^2 - X3) - S1*H^3
	let r_y = sub(mul(r, sub(dd, r_x)), mul(s1, cubed(h)));

	// Z3 = H*Z1*Z2
	let r_z = mul(mul(h, a_z), b_z);

	(r_x, r_y, r_z)
}

#[no_mangle]
pub unsafe extern "C" fn GFp_nistz256_point_add_affine(
    r: *mut Limb/*[p256::COMMON_OPS.num_limbs*3]*/,
    a: *const Limb/*[p256::COMMON_OPS.num_limbs*3]*/,
    b: *const Limb/*[p256::COMMON_OPS.num_limbs*2]*/)
{
    let a = a as *const [U256; 3];
    let b = b as *const [U256; 2];
    let r = r as *mut [U256; 3];

	// TODO: use constant
	let r_inv = U256::from_dec_str("115792089183396302114378112356516095823261736990586219612555396166510339686400").unwrap();

	let b_x: U256 = (*b)[0];
	let b_y: U256 = (*b)[1];

	if b_x.is_zero() && b_y.is_zero() {
		// a + inf = a
		*r = *a;
		return;
	}

	// b_x = b_x * r_inv^2 (mod N)
	let b_x: U256 = (((b_x.full_mul(r_inv) % U512::from(N)) * U512::from(r_inv)) % U512::from(N)).into();
	// b_x = b_x * r_inv^3 (mod N)
	let b_y: U256 = ((((b_y.full_mul(r_inv) % U512::from(N)) * U512::from(r_inv)) % U512::from(N)) * U512::from(r_inv) % U512::from(N))
		.into();

	if (*a)[2].is_zero() {
		// inf + b = b (in jacobian)
		*r = [b_x, b_y, U256::one()];
		return;
	}

    let (r_x, r_y, r_z) = p_add((*a)[0], (*a)[1], (*a)[2], b_x, b_y, U256::from(1));

    *r = [r_x, r_y, r_z];
}

#[no_mangle]
pub unsafe extern "C" fn GFp_nistz256_point_double(r: *mut Limb/*[p256::COMMON_OPS.num_limbs*3]*/,
                                                   a: *const Limb/*[p256::COMMON_OPS.num_limbs*3]*/)
{
    GFp_nistz256_point_add(r, a, a);
}
