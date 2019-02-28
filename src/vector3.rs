// vector3.rs --- Implementation of the Vector3 class, a three-dimensional vector that implements a
// number of mathematical operations.
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

use std::ops::{Add, Mul, Neg, Sub};

/// Allow a Vector3 to contain any reasonable float type.
pub trait Vector3Float:
    Default
    + Copy
    + PartialOrd
    + Add<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + Sub<Output = Self>
{
}

impl<T> Vector3Float for T where
    T: Default
        + Copy
        + PartialOrd
        + Add<Output = Self>
        + Mul<Output = Self>
        + Neg<Output = Self>
        + Sub<Output = Self>
{
}

/// A three-dimensional vector that implements a number of mathematical operations.
#[derive(Debug, Default)]
pub struct Vector3<Real: Vector3Float> {
    pub x: Real,
    pub y: Real,
    pub z: Real,
}

impl<Real: Vector3Float> Vector3<Real> {
    /// Get the dot product of two vectors.
    // TODO constexpr
    fn dot(a: &Vector3<Real>, b: &Vector3<Real>) -> Real {
        a.x * b.x + a.y * b.y + a.z * b.z
    }
}

/// The location of a plane relative to a polyhedron.
#[derive(Debug)]
pub enum PlaneLocation {
    Outside,
    Incident,
    Inside,
}

/// A two-dimensional plane defined by an offset and a normal vector.
#[derive(Debug, Default)]
pub struct Plane<Real: Vector3Float> {
    pub unit_normal: Vector3<Real>,
    pub plane_offset: Real,
}

impl<Real: Vector3Float> Plane<Real> {
    // TODO: This can be constexpr.
    fn compute_location(signed_distance: Real, tolerance: Real) -> PlaneLocation {
        if signed_distance > tolerance {
            PlaneLocation::Outside
        } else if signed_distance < -tolerance {
            PlaneLocation::Inside
        } else {
            PlaneLocation::Incident
        }
    }

    // TODO constexpr
    fn location(&self, vector: Vector3<Real>, tolerance: Real) -> PlaneLocation {
        Self::compute_location(self.signed_distance(vector), tolerance)
    }

    // TODO constexpr
    fn signed_distance(&self, vector: Vector3<Real>) -> Real {
        self.offset(vector) - self.plane_offset
    }

    // TODO constexpr
    fn offset(&self, vector: Vector3<Real>) -> Real {
        Vector3::dot(&self.unit_normal, &vector)
    }
}

#[cfg(test)]
mod tests {
    use super::Vector3;

    #[test]
    fn create_vector3() {
        let _vector = Vector3::<f64>::default();
    }
}
