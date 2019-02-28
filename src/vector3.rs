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
    // TODO constexpr
    /// Get the dot product of two vectors.
    fn dot(a: &Vector3<Real>, b: &Vector3<Real>) -> Real {
        a.x * b.x + a.y * b.y + a.z * b.z
    }
}

/// The location of a point relative to this plane.
#[derive(Debug, PartialEq)]
pub enum PlaneLocation {
    /// The point lies in the direction of the unit normal vector relative to the plane.
    Outside,
    /// The point lies on the plane.
    Incident,
    /// The point lies in the opposite direction of the unit normal vector relative to the plane.
    Inside,
}

/// A two-dimensional plane defined by an offset and a normal vector. The plane is given by
///   unit_normal.x * x + unit_normal.y * y + unit_normal.z * z + plane_offset = 0
#[derive(Debug, Default)]
pub struct Plane<Real: Vector3Float> {
    /// A unit normal vector to the plane. Note that this gives the "direction" of the plane.
    pub unit_normal: Vector3<Real>,

    /// The offset of the plane relative to the coordinate origin.
    pub plane_offset: Real,
}

impl<Real: Vector3Float> Plane<Real> {
    // TODO: This can be constexpr.
    /// Get the location of a point relative to the plane given the signed distance between the
    /// two (point - plane).
    fn location(signed_distance: Real, tolerance: Real) -> PlaneLocation {
        if signed_distance > tolerance {
            PlaneLocation::Outside
        } else if signed_distance < -tolerance {
            PlaneLocation::Inside
        } else {
            PlaneLocation::Incident
        }
    }

    // TODO constexpr
    /// Get the location of a vector relative to the plane. See the PlaneLocation struct for more
    /// details about what this means.
    pub fn vector_location(&self, vector: &Vector3<Real>, tolerance: Real) -> PlaneLocation {
        Self::location(self.signed_distance(vector), tolerance)
    }

    // TODO constexpr
    /// Get the shortest distance between a point on the plane and a vector. In this case, we
    /// simply compute the offset of a plane with the same orientation that contains the vector,
    /// and find the difference between the offsets.
    fn signed_distance(&self, vector: &Vector3<Real>) -> Real {
        self.offset_inverse(vector) - self.plane_offset
    }

    // TODO constexpr
    /// Get the the plane offset for the plane with the same unit normal vector that contains the
    /// point given by the given vector.
    ///
    /// Easy formula proof: Let (a, b, c) be the unit normal vector, and (d, e, f) be the given
    /// vector. Since (a, b, c) is a unit vector, a plane with the same unit normal vector can be
    /// given by the equation ax + by + cz + D = 0 for any distance D. Our goal is to create such a
    /// plane that contains (d, e, f). To do this, we can enter the point and solve for D.
    ///
    ///   ad + be + cf - D = 0
    ///   ad + be + cf = D
    ///
    /// Thus D = (a, b, c) * (d, e, f).
    fn offset_inverse(&self, vector: &Vector3<Real>) -> Real {
        Vector3::dot(&self.unit_normal, vector)
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
