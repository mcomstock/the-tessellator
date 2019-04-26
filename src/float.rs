// float.rs --- Implementation of the Float trait, which allows multiple types of floats to be used
// in a diagram.
// Copyright (C) 2018-2019 Maxfield Comstock

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

use crate::celery::ToCeleryPoint;
use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Neg, Sub};
use std::fmt::Debug;

pub trait CubeRoot {
    fn cbrt(self) -> Self;
}

pub trait Float:
    Default
    + Copy
    + Clone
    + Ord
    + PartialOrd
    + Add<Output = Self>
    + Div<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + Sub<Output = Self>
    + From<usize>
    + Into<usize>
    + From<f64>
    + From<i32>
    + CubeRoot
    + Debug
{
}

impl<T> Float for T where
    T: Default
        + Copy
        + Clone
        + Ord
        + PartialOrd
        + Add<Output = Self>
        + Div<Output = Self>
        + Mul<Output = Self>
        + Neg<Output = Self>
        + Sub<Output = Self>
        + From<usize>
        + Into<usize>
        + From<f64>
        + From<i32>
        + CubeRoot
        + Debug
{
}

pub trait Particle<FloatType>: ToCeleryPoint<FloatType> + Default + Clone {}

impl<T, F> Particle<F> for T where T: ToCeleryPoint<F> + Default + Clone {}

#[derive(Debug, Default, PartialOrd, PartialEq, Copy, Clone)]
pub struct Float64(pub f64);

impl Add for Float64 {
    type Output = Float64;

    fn add(self, other: Float64) -> Float64 {
        Float64(self.0 + other.0)
    }
}

impl Sub for Float64 {
    type Output = Float64;

    fn sub(self, other: Float64) -> Float64 {
        Float64(self.0 - other.0)
    }
}

impl Mul for Float64 {
    type Output = Float64;

    fn mul(self, other: Float64) -> Float64 {
        Float64(self.0 * other.0)
    }
}

impl Div for Float64 {
    type Output = Float64;

    fn div(self, other: Float64) -> Float64 {
        Float64(self.0 / other.0)
    }
}

impl Neg for Float64 {
    type Output = Float64;

    fn neg(self) -> Float64 {
        Float64(-self.0)
    }
}

impl CubeRoot for Float64 {
    fn cbrt(self) -> Float64 {
        Float64(self.0.cbrt())
    }
}

impl From<usize> for Float64 {
    fn from(u: usize) -> Self {
        Float64(u as f64)
    }
}

impl Into<usize> for Float64 {
    fn into(self) -> usize {
        self.0 as usize
    }
}

impl From<f64> for Float64 {
    fn from(f: f64) -> Self {
        Float64(f)
    }
}

impl From<i32> for Float64 {
    fn from(i: i32) -> Self {
        Float64(i as f64)
    }
}

impl Eq for Float64 {}

impl Ord for Float64 {
    fn cmp(&self, other: &Float64) -> Ordering {
        // TODO: I completely flaunt the IEEE standard for floats here, but I need a total ordering
        // to do the search. I might want to revisit this logic, or somehow ensure that having a
        // value of NaN is not possible. It seems this can be done with a wrapper type for f64.
        if self.0.is_nan() {
            Ordering::Less
        } else if other.0.is_nan() {
            Ordering::Greater
        } else {
            self.0.partial_cmp(&other.0).unwrap()
        }
    }
}
