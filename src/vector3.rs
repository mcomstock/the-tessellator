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

use crate::celery::ToCeleryPoint;
use crate::float::Float;
use std::ops::{Add, Sub};

/// A three-dimensional vector that implements a number of mathematical operations.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Vector3<Real: Float> {
    pub x: Real,
    pub y: Real,
    pub z: Real,
}

impl<Real: Float> Vector3<Real> {
    // TODO constexpr
    /// Get the dot product of two vectors.
    pub fn dot(a: &Vector3<Real>, b: &Vector3<Real>) -> Real {
        a.x * b.x + a.y * b.y + a.z * b.z
    }

    /// Get the result of scaling a vector.
    pub fn scale(&self, scalar: Real) -> Vector3<Real> {
        Vector3 {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }
}

impl<Real: Float> Add for &Vector3<Real> {
    type Output = Vector3<Real>;

    fn add(self, other: &Vector3<Real>) -> Vector3<Real> {
        Vector3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl<Real: Float> Sub for &Vector3<Real> {
    type Output = Vector3<Real>;

    fn sub(self, other: &Vector3<Real>) -> Vector3<Real> {
        Vector3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl<Real: Float> ToCeleryPoint<Real> for Vector3<Real> {
    fn get_x(&self) -> Real {
        self.x
    }

    fn get_y(&self) -> Real {
        self.y
    }

    fn get_z(&self) -> Real {
        self.z
    }
}

/// The location of a point relative to this plane.
#[derive(Clone, Debug, PartialEq)]
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
pub struct Plane<Real: Float> {
    /// A unit normal vector to the plane. Note that this gives the "direction" of the plane.
    pub unit_normal: Vector3<Real>,

    /// The offset of the plane relative to the coordinate origin.
    pub plane_offset: Real,
}

impl<Real: Float> Plane<Real> {
    // TODO: This is constexpr in the C++ version
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

    /// Get the point at the intesection of the plane and the line that passes between two points.
    pub fn intersection(&self, a: &Vector3<Real>, b: &Vector3<Real>) -> Vector3<Real> {
        let a_offset = self.offset_inverse(a);
        let b_offset = self.offset_inverse(b);

        // TODO: Make sure that this is the best way to calculate the intersection.
        return a + &(b - a).scale((self.plane_offset - a_offset) / (b_offset - a_offset));
    }
}

/// A bounding box for a space.
#[derive(Debug, Default)]
pub struct BoundingBox<Real: Float> {
    /// The low coordinates of the bounding box.
    pub low: Vector3<Real>,
    /// The high coordinates of the bounding box.
    pub high: Vector3<Real>,
}

impl<Real: Float> BoundingBox<Real> {
    /// Adjust the bounding box to contain a point.
    pub fn adjust_to_contain(&mut self, x: Real, y: Real, z: Real) {
        if x < self.low.x {
            self.low.x = x
        };

        if y < self.low.y {
            self.low.y = y
        };

        if z < self.low.z {
            self.low.z = z
        };

        if x > self.high.x {
            self.high.x = x
        };

        if y > self.high.y {
            self.high.y = y
        };

        if z > self.high.z {
            self.high.z = z
        };
    }

    /// Add padding to the bounding box.
    pub fn pad(&mut self, padding: Real) {
        self.low.x = self.low.x - padding;
        self.low.y = self.low.y - padding;
        self.low.z = self.low.z - padding;

        self.high.x = self.high.x + padding;
        self.high.y = self.high.y + padding;
        self.high.z = self.high.z + padding;
    }
}

#[cfg(test)]
mod tests {
    use super::ToCeleryPoint;
    use super::{BoundingBox, Plane, PlaneLocation, Vector3};
    use crate::float::Float64;

    fn v(x: f64, y: f64, z: f64) -> Vector3<Float64> {
        Vector3 {
            x: Float64(x),
            y: Float64(y),
            z: Float64(z),
        }
    }

    #[test]
    fn dot_test() {
        let v1 = v(1.0, 2.0, 3.0);
        let v2 = v(4.0, 5.0, 6.0);

        assert_eq!(Vector3::dot(&v1, &v2), Float64(32.0));
        assert_eq!(Vector3::dot(&v2, &v1), Float64(32.0));
    }

    #[test]
    fn scale_test() {
        let v1 = v(1.0, -2.0, 0.0);

        assert_eq!(v1.scale(Float64(1.0)), v(1.0, -2.0, 0.0));
        assert_eq!(v1.scale(Float64(2.0)), v(2.0, -4.0, 0.0));
        assert_eq!(v1.scale(Float64(-3.0)), v(-3.0, 6.0, 0.0));
    }

    #[test]
    fn add_test() {
        let v1 = v(-4.5, 0.0, 200.1);
        let v2 = v(2.0, -3.4, 4.1);

        assert_eq!(&v1 + &v2, v(-2.5, -3.4, 204.2));
        assert_eq!(&v2 + &v1, v(-2.5, -3.4, 204.2));
    }

    #[test]
    fn sub_test() {
        let v1 = v(-4.5, 0.0, 200.1);
        let v2 = v(2.0, -3.4, 4.1);

        assert_eq!(&v1 - &v2, v(-6.5, 3.4, 196.0));
        assert_eq!(&v2 - &v1, v(6.5, -3.4, -196.0));
    }

    #[test]
    fn location_test() {
        let tol = Float64(0.01);
        let out_distance = Float64(1.0);
        let in_distance = Float64(-1.0);
        let within_tolerance = Float64(0.005);

        assert_eq!(Plane::location(out_distance, tol), PlaneLocation::Outside);
        assert_eq!(Plane::location(in_distance, tol), PlaneLocation::Inside);
        assert_eq!(
            Plane::location(within_tolerance, tol),
            PlaneLocation::Incident
        );
    }

    #[test]
    fn vector_location_test() {
        // Plane that lies on x=1.
        let p = Plane {
            unit_normal: v(1.0, 0.0, 0.0),
            plane_offset: Float64(1.0),
        };

        let tol = Float64(0.05);

        let outside = v(4.0, 2.0, -6.0);
        let inside = v(-3.0, -3.0, -3.0);
        let incident = v(1.0, -6.0, 0.0);

        assert_eq!(p.vector_location(&outside, tol), PlaneLocation::Outside);
        assert_eq!(p.vector_location(&inside, tol), PlaneLocation::Inside);
        assert_eq!(p.vector_location(&incident, tol), PlaneLocation::Incident);
    }

    #[test]
    fn intersection_test() {
        // Plane that lies on x=1.
        let p = Plane {
            unit_normal: v(1.0, 0.0, 0.0),
            plane_offset: Float64(1.0),
        };

        let v1 = v(20.0, 0.0, 0.0);
        let v2 = v(10.0, 10.0, 0.0);

        assert_eq!(p.intersection(&v1, &v2), v(1.0, 19.0, 0.0));
    }

    #[test]
    fn adjust_to_contain_test() {
        let mut b = BoundingBox {
            low: v(0.0, 0.0, 0.0),
            high: v(0.0, 0.0, 0.0),
        };

        let v1 = v(-1.0, 2.0, 0.0);
        let v2 = v(-2.0, -3.0, 1.0);
        let v3 = v(0.0, 0.0, 0.0);

        b.adjust_to_contain(v1.get_x(), v1.get_y(), v1.get_z());

        assert_eq!(b.low, v(-1.0, 0.0, 0.0));
        assert_eq!(b.high, v(0.0, 2.0, 0.0));

        b.adjust_to_contain(v2.get_x(), v2.get_y(), v2.get_z());

        assert_eq!(b.low, v(-2.0, -3.0, 0.0));
        assert_eq!(b.high, v(0.0, 2.0, 1.0));

        b.adjust_to_contain(v3.get_x(), v3.get_y(), v3.get_z());

        assert_eq!(b.low, v(-2.0, -3.0, 0.0));
        assert_eq!(b.high, v(0.0, 2.0, 1.0));
    }

    #[test]
    fn pad_test() {
        let mut b = BoundingBox {
            low: v(0.0, 0.0, 0.0),
            high: v(0.0, 0.0, 0.0),
        };

        b.pad(Float64(0.5));

        assert_eq!(b.low, v(-0.5, -0.5, -0.5));
        assert_eq!(b.high, v(0.5, 0.5, 0.5));
    }
}
