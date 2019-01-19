// polyhedron.rs --- Implementation of the Polyhedron class, which performs the actual Voronoi
// calculations.
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

use pool::Pool;
use self::{StartingEdges::*, StartingFaces::*, StartingVertices::*};
use std::slice::Iter;
use vector3::Vector3;

type Vector = Vector3<f64>;

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
    point_index: Option<usize>,

    /// The index of the first half-edge in the face.
    starting_edge_index: usize,
}

/// Contains a face and additional information about that face.
#[derive(Debug, Default)]
struct FaceData {
    /// The index of the relevant face.
    face_index: usize,

    /// The vector normal to the face.
    weighted_normal: Vector,
}

// In general, the faces will be referred to using single characters,
// which are:
//   Front         : F
//   Right         : R
//   Back          : B
//   Left          : L
//   Up (Top)      : U
//   Down (Bottom) : D
#[derive(Copy, Clone, Debug)]
enum StartingFaces {
    F = 0,
    R = 1,
    B = 2,
    L = 3,
    U = 4,
    D = 5,
}

// Each vertex is written as the three faces it touches.
// The overall layout of the cube is:
//
//    BUL-------BUR               +z
//    /|        /|                |   +y
//   / |       / |                |  /
// FUL-------FUR |                | /
//  |  |      |  |       -x_______|/_______+x
//  | BDL-----|-BDR              /|
//  | /       | /               / |
//  |/        |/               /  |
// FDL-------FDR             -y   |
//                                -z
#[derive(Copy, Clone, Debug)]
enum StartingVertices {
    FDL = 0,
    FDR = 1,
    FUR = 2,
    FUL = 3,
    BDL = 4,
    BDR = 5,
    BUR = 6,
    BUL = 7,
}

#[derive(Copy, Clone, Debug)]
enum StartingEdges {
    // Front Face
    //            FU
    //    FUL <--------- FUR
    //     |              ^
    //     |              |
    //  FL |              | FR
    //     |              |
    //     v              |
    //    FDL ---------> FDR
    //            FD
    FU = 0,
    FL = 1,
    FD = 2,
    FR = 3,

    // Right Face
    //            RU
    //    FUR <--------- BUR
    //     |              ^
    //     |              |
    //  RF |              | RB
    //     |              |
    //     v              |
    //    FDR ---------> BDR
    //            RD
    RU = 4,
    RF = 5,
    RD = 6,
    RB = 7,

    // Back Face
    //            BU
    //    BUR <--------- BUL
    //     |              ^
    //     |              |
    //  BR |              | BL
    //     |              |
    //     v              |
    //    BDR ---------> BDL
    //            BD
    BU = 8,
    BR = 9,
    BD = 10,
    BL = 11,

    // Left Face
    //            LU
    //    BUL <--------- FUL
    //     |              ^
    //     |              |
    //  LB |              | LF
    //     |              |
    //     v              |
    //    BDL ---------> FDL
    //            LD
    LU = 12,
    LB = 13,
    LD = 14,
    LF = 15,

    // Up Face
    //            UB
    //    BUL <--------- BUR
    //     |              ^
    //     |              |
    //  UL |              | UR
    //     |              |
    //     v              |
    //    FUL ---------> FUR
    //            UF
    UF = 16,
    UR = 17,
    UB = 18,
    UL = 19,

    // Down Face
    //            DB
    //    BDL ---------> BDR
    //     ^              |
    //     |              |
    //  DL |              | DR
    //     |              |
    //     |              V
    //    FDL <--------- FDR
    //            DF
    DF = 20,
    DL = 21,
    DB = 22,
    DR = 23,
}

/// A cell created as the result of Voronoi tessellation. Most of the computation for the
/// tessellation is done by this struct.
#[derive(Debug, Default)]
struct Polyhedron {
    /// The index of the root edge of the polyhedron.
    root_edge: Option<usize>,

    /// The edges of the polyhedron.
    edges: Pool<HalfEdge>,
    /// The vertices around the polyhedron.
    vertices: Pool<Vector>,
    /// The faces of the polyhedron.
    faces: Pool<Face>,

    /// The face data.
    face_data: Vec<FaceData>,

    /// A list of vertices to be deleted.
    vertices_to_destroy: Vec<usize>,
    /// A list of edges do be deleted.
    edges_to_destroy: Vec<usize>,

    /// The farthest vertex from the point.
    max_distance_vertex: Option<usize>,

    /// The farthest another point can be from the current point and still cut the polyhedron.
    max_neighbor_distance: f64,
}

impl Polyhedron {
    /// Construct the initial polyhedron as a cube.
    fn build_cube(
        &mut self,
        x_min: f64,
        y_min: f64,
        z_min: f64,
        x_max: f64,
        y_max: f64,
        z_max: f64,
    ) {
        let mut vertices = Pool::<Vector>::with_capacity(8);
        let mut faces = Pool::<Face>::with_capacity(6);
        let mut edges = Pool::<HalfEdge>::with_capacity(24);

        let make_vertex = |x, y, z| Vector { x, y, z };

        // Add each vertex in order, so that the index matches the one assigned in StartingVertices.
        let fdl = vertices.add(make_vertex(x_min, y_min, z_min));
        let fdr = vertices.add(make_vertex(x_max, y_min, z_min));
        let fur = vertices.add(make_vertex(x_max, y_min, z_max));
        let ful = vertices.add(make_vertex(x_min, y_min, z_max));
        let bdl = vertices.add(make_vertex(x_min, y_max, z_min));
        let bdr = vertices.add(make_vertex(x_max, y_max, z_min));
        let bur = vertices.add(make_vertex(x_max, y_max, z_max));
        let bul = vertices.add(make_vertex(x_min, y_max, z_max));

        debug_assert_eq!(fdl, FDL as usize);
        debug_assert_eq!(fdr, FDR as usize);
        debug_assert_eq!(fur, FUR as usize);
        debug_assert_eq!(ful, FUL as usize);
        debug_assert_eq!(bdl, BDL as usize);
        debug_assert_eq!(bdr, BDR as usize);
        debug_assert_eq!(bur, BUR as usize);
        debug_assert_eq!(bul, BUL as usize);

        let make_face = |starting_edge_index| Face {
            point_index: Option::None,
            starting_edge_index: starting_edge_index as usize,
        };

        let make_edge = |face, flip, target, next| HalfEdge {
            flip: Option::Some(flip as usize),
            next: Option::Some(next as usize),
            target: Option::Some(target as usize),
            face: Option::Some(face as usize),
        };

        // Add each face in order, so that the index matches the one assigned in StartingFaces.
        let f = faces.add(make_face(FU));
        let fu = edges.add(make_edge(F, UR, FUR, RF));
        // TODO add more edges here.

        debug_assert_eq!(f, F as usize);
        debug_assert_eq!(fu, FU as usize);

        let r = faces.add(make_face(RU));

        debug_assert_eq!(r, R as usize);

        let b = faces.add(make_face(BU));

        debug_assert_eq!(b, B as usize);

        let l = faces.add(make_face(LU));

        debug_assert_eq!(l, L as usize);

        let u = faces.add(make_face(UF));

        debug_assert_eq!(u, U as usize);

        let d = faces.add(make_face(DF));

        debug_assert_eq!(d, D as usize);

        self.vertices = vertices;
        self.faces = faces;
        self.edges = edges;
    }
}

#[cfg(test)]
mod tests {
    use super::HalfEdge;

    #[test]
    fn create_half_edge() {
        let _half_edge = HalfEdge::default();
    }
}
