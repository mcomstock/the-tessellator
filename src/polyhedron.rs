// polyhedron.rs --- Implementation of the polyhedron class, which performs the actual Voronoi
// calculations
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

#[cfg(test)]
mod tests {
    use super::HalfEdge;

    #[test]
    fn create_half_edge() {
        let _half_edge = HalfEdge::default();
    }
}
