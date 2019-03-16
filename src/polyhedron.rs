// polyhedron.rs --- Implementation of the polyhedron class, which performs the actual Voronoi
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

use self::{StartingEdges::*, StartingFaces::*, StartingVertices::*};
use crate::pool::Pool;
use crate::vector3::{Plane, PlaneLocation, Vector3};

use std::cell::RefCell;

type Vector = Vector3<f64>;
type Plane64 = Plane<f64>;

const TOLERANCE: f64 = 1e-12;

/// A half-edge structure. Two half-edges make up an edge in the polyhedron, each pointing to one
/// of the two vertices that make up that edge.
#[derive(Debug, Default)]
struct HalfEdge {
    /// The index of the other half-edge that makes up this edge.
    flip: Option<usize>,

    /// The index of the next half-edge in the polyhedron face.
    next: Option<usize>,

    // TODO: This might not need interior mutability to work.
    /// The vertex this half-edge belongs to.
    target: RefCell<Option<usize>>,

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
    /// Get a new polyhedron initialized as a cube.
    pub fn new(
        x_min: f64,
        y_min: f64,
        z_min: f64,
        x_max: f64,
        y_max: f64,
        z_max: f64,
    ) -> Polyhedron {
        let (root_edge, vertices, faces, edges) =
            Polyhedron::build_cube(x_min, y_min, z_min, x_max, y_max, z_max);

        Polyhedron {
            root_edge,
            vertices,
            faces,
            edges,
            ..Default::default()
        }
    }

    /// Construct the initial polyhedron as a cube.
    fn build_cube(
        x_min: f64,
        y_min: f64,
        z_min: f64,
        x_max: f64,
        y_max: f64,
        z_max: f64,
    ) -> (Option<usize>, Pool<Vector>, Pool<Face>, Pool<HalfEdge>) {
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
            point_index: None,
            starting_edge_index: starting_edge_index as usize,
        };

        let make_edge = |face, flip, target, next| HalfEdge {
            flip: Some(flip as usize),
            next: Some(next as usize),
            target: RefCell::new(Some(target as usize)),
            face: Some(face as usize),
        };

        // Add each face in order, so that the index matches the one assigned in StartingFaces.
        let f = faces.add(make_face(FU));
        let fu = edges.add(make_edge(F, UF, FUL, FL));
        let fl = edges.add(make_edge(F, LF, FDL, FD));
        let fd = edges.add(make_edge(F, DF, FDR, FR));
        let fr = edges.add(make_edge(F, RF, FUR, FU));

        debug_assert_eq!(f, F as usize);
        debug_assert_eq!(fu, FU as usize);
        debug_assert_eq!(fl, FL as usize);
        debug_assert_eq!(fd, FD as usize);
        debug_assert_eq!(fr, FR as usize);

        let r = faces.add(make_face(RU));
        let ru = edges.add(make_edge(R, UR, FUR, RF));
        let rf = edges.add(make_edge(R, FR, FDR, RD));
        let rd = edges.add(make_edge(R, DR, BDR, RB));
        let rb = edges.add(make_edge(R, BR, BUR, RU));

        debug_assert_eq!(r, R as usize);
        debug_assert_eq!(ru, RU as usize);
        debug_assert_eq!(rf, RF as usize);
        debug_assert_eq!(rd, RD as usize);
        debug_assert_eq!(rb, RB as usize);

        let b = faces.add(make_face(BU));
        let bu = edges.add(make_edge(B, UB, BUR, BR));
        let br = edges.add(make_edge(B, RB, BDR, BD));
        let bd = edges.add(make_edge(B, DB, BDL, BL));
        let bl = edges.add(make_edge(B, LB, BUL, BU));

        debug_assert_eq!(b, B as usize);
        debug_assert_eq!(bu, BU as usize);
        debug_assert_eq!(br, BR as usize);
        debug_assert_eq!(bd, BD as usize);
        debug_assert_eq!(bl, BL as usize);

        let l = faces.add(make_face(LU));
        let lu = edges.add(make_edge(L, UL, BUL, LB));
        let lb = edges.add(make_edge(L, BL, BDL, LD));
        let ld = edges.add(make_edge(L, DL, FDL, LF));
        let lf = edges.add(make_edge(L, FL, FUL, LU));

        debug_assert_eq!(l, L as usize);
        debug_assert_eq!(lu, LU as usize);
        debug_assert_eq!(lb, LB as usize);
        debug_assert_eq!(ld, LD as usize);
        debug_assert_eq!(lf, LF as usize);

        let u = faces.add(make_face(UF));
        let uf = edges.add(make_edge(U, FU, FUR, UR));
        let ur = edges.add(make_edge(U, RU, BUR, UB));
        let ub = edges.add(make_edge(U, BU, BUL, UL));
        let ul = edges.add(make_edge(U, LU, FUL, UF));

        debug_assert_eq!(u, U as usize);
        debug_assert_eq!(uf, UF as usize);
        debug_assert_eq!(ur, UR as usize);
        debug_assert_eq!(ub, UB as usize);
        debug_assert_eq!(ul, UL as usize);

        let d = faces.add(make_face(DF));
        let df = edges.add(make_edge(D, FD, FDL, DL));
        let dl = edges.add(make_edge(D, LD, BDL, DB));
        let db = edges.add(make_edge(D, BD, BDR, DR));
        let dr = edges.add(make_edge(B, RD, FDR, DF));

        debug_assert_eq!(d, D as usize);
        debug_assert_eq!(df, DF as usize);
        debug_assert_eq!(dl, DL as usize);
        debug_assert_eq!(db, DB as usize);
        debug_assert_eq!(dr, DR as usize);

        (Some(FU as usize), vertices, faces, edges)
    }

    /// Find the index of an edge that is going out of the plane. The edge does not need to begin
    /// inside the plane to satisfy this condition.
    fn find_outgoing_edge(&self, plane: &Plane64) -> Option<usize> {
        // Check for any vertices that would be cut off by the plane - these are "outside" vertices.
        // If no such vertices are present, then there is no need to cut at all.
        let mut need_to_cut = false;
        for vertex in self.vertices.into_iter() {
            if plane.vector_location(vertex, TOLERANCE) == PlaneLocation::Outside {
                need_to_cut = true;
                break;
            }
        }

        // If no outside edges were found, simply return nothing.
        if !need_to_cut {
            return None;
        }

        // Search the edges until an outgoing one is found.
        for edge in self.edges.into_iter() {
            let target_index = edge.target.borrow().unwrap();
            let target = self.vertices.get(target_index).unwrap();

            // If the edge is ingoing, check its flip.
            // TODO check if this is actually faster
            match plane.vector_location(&target, TOLERANCE) {
                PlaneLocation::Inside => {
                    let flip_target_index = self.target_index(edge.flip).unwrap();
                    let flip_target = self.vertices.get(flip_target_index).unwrap();

                    if plane.vector_location(flip_target, TOLERANCE) == PlaneLocation::Outside {
                        return Some(flip_target_index);
                    }
                }

                // TODO the C++ version doesn't check this case - why?
                PlaneLocation::Outside => {
                    return Some(target_index);
                }

                _ => continue,
            }
        }

        return None;
    }

    /// Cut the face of the polyhedron with the given plane.
    fn cut_with_plane(&mut self, face_to_cut_index: usize, plane: &Plane64) -> bool {
        // TODO check validity for debug?

        let first_outgoing_edge_index = match self.find_outgoing_edge(plane) {
            Some(index) => {
                // Check that the target vertex is not inside.
                debug_assert_ne!(
                    plane.vector_location(
                        &self
                            .vertices
                            .get_or_fail(self.edges.get_or_fail(index).target.borrow().unwrap()),
                        TOLERANCE
                    ),
                    PlaneLocation::Inside
                );

                // Check that the source vertex (target of the flip) is inside.
                debug_assert_eq!(
                    plane.vector_location(
                        &self.vertices.get_or_fail(
                            self.target_index(self.edges.get_or_fail(index).flip)
                                .unwrap()
                        ),
                        TOLERANCE
                    ),
                    PlaneLocation::Inside
                );

                index
            }
            // This is a legitimate case where there is no outgoing edge to cut.
            _ => return false,
        };

        let first_outside_face_edge_index = self.edges.next_index();
        let outside_face_index = self.faces.next_index();

        let actual_edge_index = self.edges.add(HalfEdge {
            face: Some(outside_face_index),
            ..Default::default()
        });

        let actual_face_index = self.faces.add(Face {
            // TODO this doesn't seem right
            point_index: Some(face_to_cut_index),
            starting_edge_index: first_outside_face_edge_index,
        });

        debug_assert_eq!(first_outside_face_edge_index, actual_edge_index);
        debug_assert_eq!(outside_face_index, actual_face_index);

        // TODO check if it's fine to do this on-demand instead of batching
        let mut vertices_to_destroy = Vec::<usize>::new();
        let mut edges_to_destroy = Vec::<usize>::new();

        let mut outside_face_edge_index = first_outside_face_edge_index;
        let mut outgoing_edge_index = first_outgoing_edge_index;

        let mut previous_intersection: Option<usize> = None;

        loop {
            let mut previous_vertex_index =
                self.target_index_from_index(outgoing_edge_index).unwrap();
            vertices_to_destroy.push(previous_vertex_index);

            // The target vertex should not be inside.
            debug_assert_ne!(
                plane.vector_location(self.vertices.get_or_fail(previous_vertex_index), TOLERANCE),
                PlaneLocation::Inside
            );

            // Set the next edge and target vertex to potentially cut.
            let mut current_edge_index = self.next_index_from_index(outgoing_edge_index).unwrap();
            let mut current_vertex_index =
                self.target_index_from_index(current_edge_index).unwrap();

            // TODO: Is there a better way to check this?
            // If the starting vertex was outside the plane, we definitely need to cut.
            let mut need_to_cut = plane
                .vector_location(self.vertices.get(previous_vertex_index).unwrap(), TOLERANCE)
                == PlaneLocation::Outside;

            let mut previous_location =
                plane.vector_location(self.vertices.get(current_vertex_index).unwrap(), TOLERANCE);
            let mut current_location = previous_location.clone();

            // Keep going until we are back inside.
            while current_location != PlaneLocation::Inside {
                need_to_cut = true;
                vertices_to_destroy.push(current_vertex_index);
                edges_to_destroy.push(current_edge_index);

                previous_vertex_index = current_vertex_index;
                previous_location = current_location;

                current_edge_index = self.next_index_from_index(current_edge_index).unwrap();
                current_vertex_index = self.target_index_from_index(current_edge_index).unwrap();

                current_location = plane
                    .vector_location(self.vertices.get(current_vertex_index).unwrap(), TOLERANCE);
            }

            // TODO: This should cause a bug if `previous_intersection` is not set in the loop.
            self.edges
                .get_or_fail(outgoing_edge_index)
                .target
                .replace(Some(previous_intersection.unwrap()));

            if need_to_cut {
                let current_intersection_vertex_index = if previous_location == PlaneLocation::Incident {
                    let old_vertex = self.vertices.get_or_fail(previous_vertex_index).to_owned();

                    // TODO: Can this be cloned instead?
                    let vertex: Vector = Vector {
                        x: old_vertex.x,
                        y: old_vertex.y,
                        z: old_vertex.z,
                    };

                    self.vertices.add(vertex)
                } else {
                    // TODO: Build new vertex at intersection.
                    1
                };
            }

            // The test condition that should break out of the loop.
            if outgoing_edge_index == first_outgoing_edge_index {
                break;
            }
        }

        false
    }

    /// Get the target vertex index of the edge at the given edge index, if possible.
    fn target_index_from_index(&self, edge_index: usize) -> Option<usize> {
        // TODO get rid of clone
        self.edges
            .get(edge_index)
            .and_then(|e| e.target.borrow().clone())
    }

    /// Get the target vertex index of the edge at the given edge index, if possible.
    fn target_index(&self, edge_index: Option<usize>) -> Option<usize> {
        edge_index.and_then(|i| self.target_index_from_index(i))
    }

    /// Get the source vertex index of the edge at the given edge index, if possible.
    fn source_index(&self, edge_index: Option<usize>) -> Option<usize> {
        edge_index
            .and_then(|i| self.edges.get(i))
            .and_then(|e| self.target_index(e.flip))
    }

    /// Get the next edge from the current edge index, if possible.
    fn next_index_from_index(&self, edge_index: usize) -> Option<usize> {
        self.edges.get(edge_index).and_then(|e| e.next)
    }

    /// Get the next edge from the current edge index, if possible.
    fn next_index(&self, edge_index: Option<usize>) -> Option<usize> {
        edge_index.and_then(|i| self.next_index_from_index(i))
    }
}

#[cfg(test)]
mod tests {
    use super::{HalfEdge, Polyhedron};

    #[test]
    fn create_half_edge() {
        let _half_edge = HalfEdge::default();
    }

    #[test]
    fn build_cube() {
        Polyhedron::build_cube(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
    }

    #[test]
    fn new() {
        let p = Polyhedron::new(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
    }
}
