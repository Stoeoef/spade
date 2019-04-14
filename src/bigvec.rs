// Copyright 2017 The Spade Developers.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use crate::point_traits::{PointN, TwoDimensional};
use crate::traits::SpadeNum;
use cgmath as cg;
use nalgebra as na;
use num::bigint::ToBigInt;
use num::{BigInt, Integer, Num, One, Signed, ToPrimitive, Zero};
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Rem, Sub};

/// BigVec2 is a two dimensional vector that does not _require_ it's
/// internal scalar type to be `Copy`.
#[repr(C)]
#[derive(Clone, Debug, PartialEq, Eq, Copy)]
pub struct BigVec2<N: Num> {
    pub x: N,
    pub y: N,
}

impl<N: Num + Clone> BigVec2<N> {
    pub fn new(x: N, y: N) -> BigVec2<N> {
        BigVec2 { x, y }
    }
}

impl<N: SpadeNum> PointN for BigVec2<N> {
    type Scalar = N;

    fn dimensions() -> usize {
        2
    }

    fn nth(&self, index: usize) -> &N {
        &self[index]
    }
    fn nth_mut(&mut self, index: usize) -> &mut N {
        &mut self[index]
    }

    fn from_value(val: N) -> Self {
        BigVec2::new(val.clone(), val.clone())
    }
}

impl<S: SpadeNum> TwoDimensional for BigVec2<S> {}

impl<N: SpadeNum> Add for BigVec2<N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        BigVec2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<N: SpadeNum> Sub for BigVec2<N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        BigVec2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl<N: SpadeNum> Mul<N> for BigVec2<N> {
    type Output = Self;
    fn mul(self, rhs: N) -> Self {
        BigVec2::new(self.x * rhs.clone(), self.y * rhs.clone())
    }
}

impl<N: SpadeNum> Div<N> for BigVec2<N> {
    type Output = Self;
    fn div(self, rhs: N) -> Self {
        BigVec2::new(self.x / rhs.clone(), self.y / rhs.clone())
    }
}

impl<N: SpadeNum> Index<usize> for BigVec2<N> {
    type Output = N;

    fn index(&self, index: usize) -> &N {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => unreachable!(),
        }
    }
}

impl<N: SpadeNum> IndexMut<usize> for BigVec2<N> {
    fn index_mut(&mut self, index: usize) -> &mut N {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => unreachable!(),
        }
    }
}

impl<I> From<cg::Point2<I>> for BigVec2<BigInt>
where
    I: cg::BaseNum + ToBigInt,
{
    fn from(v: cg::Point2<I>) -> Self {
        BigVec2::new(v[0].to_bigint().unwrap(), v[1].to_bigint().unwrap())
    }
}

impl<I> From<na::Point2<I>> for BigVec2<BigInt>
where
    I: ToBigInt + na::Scalar,
{
    fn from(v: na::Point2<I>) -> Self {
        BigVec2::new(v[0].to_bigint().unwrap(), v[1].to_bigint().unwrap())
    }
}

#[derive(Debug, Clone)]
pub enum AdaptiveInt {
    LowRes(i64),
    HighRes(BigInt),
}

impl AdaptiveInt {
    pub fn from_i64(i: i64) -> AdaptiveInt {
        AdaptiveInt::LowRes(i)
    }

    pub fn from_bigint(i: BigInt) -> AdaptiveInt {
        AdaptiveInt::HighRes(i.clone()).reduce()
    }

    pub fn reduce(self) -> AdaptiveInt {
        match self {
            AdaptiveInt::LowRes(i) => AdaptiveInt::LowRes(i),
            AdaptiveInt::HighRes(big) => {
                if let Some(i) = big.to_i64() {
                    AdaptiveInt::LowRes(i)
                } else {
                    AdaptiveInt::HighRes(big)
                }
            }
        }
    }
}

impl Num for AdaptiveInt {
    type FromStrRadixErr = ::std::num::ParseIntError;
    fn from_str_radix(s: &str, radix: u32) -> Result<Self, ::std::num::ParseIntError> {
        i64::from_str_radix(s, radix).map(AdaptiveInt::LowRes)
    }
}

impl Integer for AdaptiveInt {
    fn div_floor(&self, other: &Self) -> Self {
        use crate::bigvec::AdaptiveInt::*;
        match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => HighRes(l.div_floor(r)).reduce(),
            (&HighRes(ref l), &LowRes(ref r)) => {
                HighRes(l.div_floor(&r.to_bigint().unwrap())).reduce()
            }
            (&LowRes(ref l), &HighRes(ref r)) => {
                HighRes(l.to_bigint().unwrap().div_floor(r)).reduce()
            }
            (&LowRes(ref l), &LowRes(ref r)) => LowRes(l.div_floor(r)),
        }
    }

    fn mod_floor(&self, other: &Self) -> Self {
        use crate::bigvec::AdaptiveInt::*;
        match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => HighRes(l.mod_floor(r)).reduce(),
            (&HighRes(ref l), &LowRes(ref r)) => {
                HighRes(l.mod_floor(&r.to_bigint().unwrap())).reduce()
            }
            (&LowRes(ref l), &HighRes(ref r)) => {
                HighRes(l.to_bigint().unwrap().mod_floor(r)).reduce()
            }
            (&LowRes(ref l), &LowRes(ref r)) => LowRes(l.mod_floor(r)),
        }
    }

    fn gcd(&self, other: &Self) -> Self {
        use crate::bigvec::AdaptiveInt::*;
        match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => HighRes(l.gcd(r)).reduce(),
            (&HighRes(ref l), &LowRes(ref r)) => HighRes(l.gcd(&r.to_bigint().unwrap())).reduce(),
            (&LowRes(ref l), &HighRes(ref r)) => HighRes(l.to_bigint().unwrap().gcd(r)).reduce(),
            (&LowRes(ref l), &LowRes(ref r)) => LowRes(l.gcd(r)),
        }
    }

    fn lcm(&self, other: &Self) -> Self {
        use crate::bigvec::AdaptiveInt::*;
        match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => HighRes(l.lcm(r)),
            (&HighRes(ref l), &LowRes(ref r)) => HighRes(l.lcm(&r.to_bigint().unwrap())),
            (&LowRes(ref l), &HighRes(ref r)) => HighRes(l.to_bigint().unwrap().lcm(r)),
            (&LowRes(ref l), &LowRes(ref r)) => LowRes(l.lcm(r)),
        }
    }

    fn divides(&self, other: &Self) -> bool {
        use crate::bigvec::AdaptiveInt::*;
        match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => l.divides(r),
            (&HighRes(ref l), &LowRes(ref r)) => l.divides(&r.to_bigint().unwrap()),
            (&LowRes(ref l), &HighRes(ref r)) => l.to_bigint().unwrap().divides(r),
            (&LowRes(ref l), &LowRes(ref r)) => l.divides(r),
        }
    }

    fn is_multiple_of(&self, other: &Self) -> bool {
        use crate::bigvec::AdaptiveInt::*;
        match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => l.is_multiple_of(r),
            (&HighRes(ref l), &LowRes(ref r)) => l.is_multiple_of(&r.to_bigint().unwrap()),
            (&LowRes(ref l), &HighRes(ref r)) => l.to_bigint().unwrap().is_multiple_of(r),
            (&LowRes(ref l), &LowRes(ref r)) => l.is_multiple_of(r),
        }
    }

    fn is_even(&self) -> bool {
        match self {
            AdaptiveInt::HighRes(ref v) => v.is_even(),
            AdaptiveInt::LowRes(ref v) => v.is_even(),
        }
    }

    fn is_odd(&self) -> bool {
        !self.is_even()
    }

    fn div_rem(&self, other: &Self) -> (Self, Self) {
        use crate::bigvec::AdaptiveInt::*;
        let (div, rem) = match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => l.div_rem(r),
            (&HighRes(ref l), &LowRes(ref r)) => l.div_rem(&r.to_bigint().unwrap()),
            (&LowRes(ref l), &HighRes(ref r)) => l.to_bigint().unwrap().div_rem(r),
            (&LowRes(ref l), &LowRes(ref r)) => {
                let (div, rem) = l.div_rem(r);
                return (LowRes(div), LowRes(rem));
            }
        };
        (HighRes(div).reduce(), HighRes(rem).reduce())
    }
}

impl Zero for AdaptiveInt {
    fn zero() -> Self {
        AdaptiveInt::LowRes(0)
    }
    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }
}

impl One for AdaptiveInt {
    fn one() -> Self {
        AdaptiveInt::LowRes(1)
    }
}

impl Add for AdaptiveInt {
    type Output = Self;
    fn add(self, rhs: AdaptiveInt) -> AdaptiveInt {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (HighRes(l), HighRes(r)) => HighRes(l + r).reduce(),
            (HighRes(hr), LowRes(lr)) | (LowRes(lr), HighRes(hr)) => {
                HighRes(hr + lr.to_bigint().unwrap()).reduce()
            }
            (LowRes(l), LowRes(r)) => {
                if let Some(sum) = l.checked_add(r) {
                    LowRes(sum)
                } else {
                    HighRes(l.to_bigint().unwrap() + r.to_bigint().unwrap())
                }
            }
        }
    }
}

impl Sub for AdaptiveInt {
    type Output = Self;
    fn sub(self, rhs: AdaptiveInt) -> AdaptiveInt {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (HighRes(l), HighRes(r)) => HighRes(l - r).reduce(),
            (HighRes(l), LowRes(r)) => HighRes(l - r.to_bigint().unwrap()).reduce(),
            (LowRes(l), HighRes(r)) => HighRes(l.to_bigint().unwrap() - r).reduce(),
            (LowRes(l), LowRes(r)) => {
                if let Some(diff) = l.checked_sub(r) {
                    LowRes(diff)
                } else {
                    HighRes(l.to_bigint().unwrap() - r.to_bigint().unwrap())
                }
            }
        }
    }
}

impl Mul for AdaptiveInt {
    type Output = Self;
    fn mul(self, rhs: AdaptiveInt) -> AdaptiveInt {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (HighRes(l), HighRes(r)) => HighRes(l * r),
            (HighRes(hr), LowRes(lr)) | (LowRes(lr), HighRes(hr)) => {
                HighRes(hr * lr.to_bigint().unwrap())
            }
            (LowRes(l), LowRes(r)) => {
                if let Some(prod) = l.checked_mul(r) {
                    LowRes(prod)
                } else {
                    HighRes(l.to_bigint().unwrap() * r.to_bigint().unwrap())
                }
            }
        }
    }
}

impl Div for AdaptiveInt {
    type Output = Self;

    fn div(self, rhs: AdaptiveInt) -> AdaptiveInt {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (HighRes(l), HighRes(r)) => HighRes(l / r).reduce(),
            (HighRes(l), LowRes(r)) => HighRes(l / r.to_bigint().unwrap()).reduce(),
            (LowRes(l), HighRes(r)) => HighRes(l.to_bigint().unwrap() / r).reduce(),
            (LowRes(l), LowRes(r)) => {
                if let Some(quot) = l.checked_div(r) {
                    LowRes(quot)
                } else {
                    HighRes(l.to_bigint().unwrap() / r.to_bigint().unwrap()).reduce()
                }
            }
        }
    }
}

impl Rem for AdaptiveInt {
    type Output = Self;
    fn rem(self, rhs: AdaptiveInt) -> AdaptiveInt {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (HighRes(l), HighRes(r)) => HighRes(l % r),
            (HighRes(l), LowRes(r)) => HighRes(l % r.to_bigint().unwrap()),
            (LowRes(l), HighRes(r)) => HighRes(l.to_bigint().unwrap() % r),
            (LowRes(l), LowRes(r)) => {
                if let Some(quot) = l.checked_rem(r) {
                    LowRes(quot)
                } else {
                    HighRes(l.to_bigint().unwrap() % r.to_bigint().unwrap())
                }
            }
        }
    }
}

impl Neg for AdaptiveInt {
    type Output = AdaptiveInt;
    fn neg(self) -> Self {
        match self {
            AdaptiveInt::HighRes(v) => AdaptiveInt::HighRes(-v),
            AdaptiveInt::LowRes(v) => AdaptiveInt::LowRes(-v),
        }
    }
}

impl Signed for AdaptiveInt {
    fn abs(&self) -> Self {
        match self {
            AdaptiveInt::HighRes(ref v) => AdaptiveInt::HighRes(v.abs()),
            AdaptiveInt::LowRes(ref v) => AdaptiveInt::LowRes(v.abs()),
        }
    }
    fn abs_sub(&self, other: &Self) -> Self {
        use crate::bigvec::AdaptiveInt::*;
        match (self, other) {
            (&HighRes(ref l), &HighRes(ref r)) => HighRes(l.abs_sub(r)),
            (&LowRes(ref l), &HighRes(ref r)) => HighRes(l.to_bigint().unwrap().abs_sub(r)),
            (&HighRes(ref l), &LowRes(ref r)) => HighRes(l.abs_sub(&r.to_bigint().unwrap())),
            (&LowRes(ref l), &LowRes(ref r)) => LowRes(l.abs_sub(r)),
        }
    }

    fn signum(&self) -> Self {
        match self {
            AdaptiveInt::HighRes(ref v) => AdaptiveInt::HighRes(v.signum()),
            AdaptiveInt::LowRes(ref v) => AdaptiveInt::LowRes(v.signum()),
        }
    }

    fn is_positive(&self) -> bool {
        match self {
            AdaptiveInt::HighRes(ref v) => v.is_positive(),
            AdaptiveInt::LowRes(ref v) => v.is_positive(),
        }
    }

    fn is_negative(&self) -> bool {
        match self {
            AdaptiveInt::HighRes(ref v) => v.is_negative(),
            AdaptiveInt::LowRes(ref v) => v.is_negative(),
        }
    }
}

impl PartialEq for AdaptiveInt {
    fn eq(&self, rhs: &AdaptiveInt) -> bool {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (&HighRes(ref l), &HighRes(ref r)) => l == r,
            (&HighRes(ref hr), &LowRes(ref lr)) | (&LowRes(ref lr), &HighRes(ref hr)) => {
                hr == &lr.to_bigint().unwrap()
            }
            (&LowRes(ref l), &LowRes(ref r)) => l == r,
        }
    }
}

impl Eq for AdaptiveInt {}

impl PartialOrd for AdaptiveInt {
    fn partial_cmp(&self, rhs: &AdaptiveInt) -> Option<::std::cmp::Ordering> {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (&HighRes(ref l), &HighRes(ref r)) => l.partial_cmp(r),
            (&HighRes(ref l), &LowRes(ref r)) => l.partial_cmp(&r.to_bigint().unwrap()),
            (&LowRes(ref l), &HighRes(ref r)) => l.to_bigint().unwrap().partial_cmp(r),
            (&LowRes(ref l), &LowRes(ref r)) => l.partial_cmp(r),
        }
    }
}

impl Ord for AdaptiveInt {
    fn cmp(&self, rhs: &AdaptiveInt) -> ::std::cmp::Ordering {
        use crate::bigvec::AdaptiveInt::*;
        match (self, rhs) {
            (&HighRes(ref l), &HighRes(ref r)) => l.cmp(r),
            (&HighRes(ref l), &LowRes(ref r)) => l.cmp(&r.to_bigint().unwrap()),
            (&LowRes(ref l), &HighRes(ref r)) => l.to_bigint().unwrap().cmp(r),
            (&LowRes(ref l), &LowRes(ref r)) => l.cmp(r),
        }
    }
}

impl<I> From<cg::Point2<I>> for BigVec2<AdaptiveInt>
where
    I: ::std::convert::Into<i64> + Copy,
{
    fn from(v: cg::Point2<I>) -> Self {
        BigVec2::new(
            AdaptiveInt::LowRes(v[0].into()),
            AdaptiveInt::LowRes(v[1].into()),
        )
    }
}

impl<I> From<na::Point2<I>> for BigVec2<AdaptiveInt>
where
    I: ::std::convert::Into<i64> + na::Scalar,
{
    fn from(v: na::Point2<I>) -> Self {
        BigVec2::new(
            AdaptiveInt::LowRes(v[0].into()),
            AdaptiveInt::LowRes(v[1].into()),
        )
    }
}
