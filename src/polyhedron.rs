// polyhedron.rs --- Implementation of the Polyhedron class, which performs the actual Voronoi
// calculations.
// Copyright (C) 2018 Maxfield Comstock

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

use vector3::Vector3;

/// A half-edge structure. Two half-edges make up an edge in the polyhedron, each pointing to one
/// of the two vertices that make up that edge.
#[derive(Debug, Default)]
struct HalfEdge {
    /// The index of the other half-edge that makes up this edge.
    flip: Option<usize>,

    /// The index of the next half-edge in the polyhedron face.
    next: Option<usize>,

    /// The vertex this half-edge belongs to.
    target: Option<usize>,

    /// The polyhedron face this half-edge is in.
    face: Option<usize>,
}

/// A face of a polyhedron.
#[derive(Debug, Default)]
struct Face {
    /// The index of the point that the face contains.
    point_index: usize,

    /// The index of the first half-edge in the face.
    starting_edge_index: usize,
}

/// Contains a face and additional information about that face.
#[derive(Debug, Default)]
struct FaceData {
    /// The index of the relevant face.
    face_index: usize,

    /// The vector normal to the face.
    weighted_normal: Vector3<f64>,
}

/// A cell created as the result of Voronoi tessellation. Most of the computation for the
/// tessellation is done by this struct.
#[derive(Debug, Default)]
struct Polyhedron {

}

#[cfg(test)]
mod tests {
    use super::HalfEdge;

    #[test]
    fn create_half_edge() {
        let _half_edge = HalfEdge::default();
    }
}
