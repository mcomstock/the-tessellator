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

/// A three-dimensional vector that implements a number of mathematical operations.
#[derive(Debug, Default)]
pub struct Vector3<Real: Default> {
    pub x: Real,
    pub y: Real,
    pub z: Real,
}

#[cfg(test)]
mod tests {
    use super::Vector3;

    #[test]
    fn create_vector3() {
        let _vector = Vector3::<f64>::default();
    }
}
