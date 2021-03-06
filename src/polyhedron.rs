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
use crate::float::Float;
use crate::pool::Pool;
use crate::vector3::{Plane, PlaneLocation, Vector3};

/// A half-edge structure. Two half-edges make up an edge in the polyhedron, each pointing to one
/// of the two vertices that make up that edge.
#[derive(Debug, Default, PartialEq)]
pub struct HalfEdge {
    /// The index of the other half-edge that makes up this edge.
    pub flip: Option<usize>,

    /// The index of the next half-edge in the polyhedron face.
    pub next: Option<usize>,

    /// The vertex this half-edge belongs to.
    pub target: Option<usize>,

    /// The polyhedron face this half-edge is in.
    pub face: Option<usize>,
}

// TODO impl and use convenience methods for getting the target, source, flip, etc

/// A face of a polyhedron.
#[derive(Debug, Default)]
pub struct Face {
    /// The index of the point that the face contains. This is not used by the Polyhedron, but
    /// is used by the external diagram structure to determine which of its points created this
    /// face.
    pub point_index: Option<usize>,

    /// The index of the first half-edge in the face.
    pub starting_edge_index: usize,
}

// TODO: Think about making this an optional field on the Face itself.
/// Contains a face and additional information about that face.
#[derive(Debug, Default)]
pub struct FaceData<Real: Float> {
    /// The index of the relevant face.
    pub face_index: usize,

    /// The vector normal to the face.
    pub weighted_normal: Vector3<Real>,
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
pub struct Polyhedron<Real: Float> {
    /// The index of the root edge of the polyhedron.
    pub root_edge: Option<usize>,

    /// The edges of the polyhedron.
    pub edges: Pool<HalfEdge>,
    /// The vertices around the polyhedron.
    pub vertices: Pool<Vector3<Real>>,
    /// The faces of the polyhedron.
    pub faces: Pool<Face>,

    /// The face data.
    pub face_data: Vec<FaceData<Real>>,
}

impl<Real: Float> Polyhedron<Real> {
    /// The tolerance within which two floats should be considered equal.
    fn tolerance() -> Real {
        Real::from(1e-12)
    }

    /// Get a new polyhedron initialized as a cube.
    pub fn new(
        x_min: Real,
        y_min: Real,
        z_min: Real,
        x_max: Real,
        y_max: Real,
        z_max: Real,
    ) -> Polyhedron<Real> {
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

    /// Reset the polyhedron to a cube
    pub fn reset(
        &mut self,
        x_min: Real,
        y_min: Real,
        z_min: Real,
        x_max: Real,
        y_max: Real,
        z_max: Real,
    ) {
        let (root_edge, vertices, faces, edges) =
            Polyhedron::build_cube(x_min, y_min, z_min, x_max, y_max, z_max);

        self.root_edge = root_edge;
        self.vertices = vertices;
        self.faces = faces;
        self.edges = edges;

        self.face_data.clear();
    }

    /// Construct the initial polyhedron as a cube.
    fn build_cube(
        x_min: Real,
        y_min: Real,
        z_min: Real,
        x_max: Real,
        y_max: Real,
        z_max: Real,
    ) -> (
        Option<usize>,
        Pool<Vector3<Real>>,
        Pool<Face>,
        Pool<HalfEdge>,
    ) {
        let mut vertices = Pool::<Vector3<Real>>::with_capacity(8);
        let mut faces = Pool::<Face>::with_capacity(6);
        let mut edges = Pool::<HalfEdge>::with_capacity(24);

        let make_vertex = |x, y, z| Vector3::<Real> { x, y, z };

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
            target: Some(target as usize),
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

    /// Find the index of an edge with its target outside of the plane and its source out of the
    /// plane.
    fn find_outgoing_edge(&self, plane: &Plane<Real>) -> Option<usize> {
        // Check for any vertices that would be cut off by the plane - these are "outside" vertices.
        // If no such vertices are present, then there is no need to cut at all.
        let mut need_to_cut = false;
        for vertex in self.vertices.into_iter() {
            if plane.vector_location(vertex, Polyhedron::tolerance()) == PlaneLocation::Outside {
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
            let target_index = edge.target.unwrap();
            let target = self.vertices.get(target_index).unwrap();

            // If the edge is ingoing, check its flip.
            match plane.vector_location(&target, Polyhedron::tolerance()) {
                PlaneLocation::Inside => {
                    let flip_target_index = self.target_index(edge.flip).unwrap();
                    let flip_target = self.vertices.get(flip_target_index).unwrap();

                    if plane.vector_location(flip_target, Polyhedron::tolerance())
                        == PlaneLocation::Outside
                    {
                        return Some(edge.flip.unwrap());
                    }
                }

                _ => continue,
            }
        }

        return None;
    }

    /// Cut the face of the polyhedron with the given plane.
    pub fn cut_with_plane(&mut self, point_index: usize, plane: &Plane<Real>) -> bool {
        // Find an edge crossing the plane, if one exists.
        let first_outgoing_edge_index = match self.find_outgoing_edge(plane) {
            Some(index) => {
                // Check that the target vertex is not inside.
                debug_assert_ne!(
                    plane.vector_location(
                        &self
                            .vertices
                            .get_or_fail(self.edges.get_or_fail(index).target.unwrap()),
                        Polyhedron::tolerance()
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
                        Polyhedron::tolerance()
                    ),
                    PlaneLocation::Inside
                );

                index
            }
            // This is a legitimate case where there is no outgoing edge to cut.
            _ => return false,
        };

        // Set the root edge to something that is certain not to be destroyed by the cutting
        // process.
        self.root_edge = Some(first_outgoing_edge_index);

        let mut outgoing_edge_index = first_outgoing_edge_index;
        let mut previous_intersection: Option<usize> = None;

        let first_outside_face_edge_index = self.edges.add(HalfEdge::default());
        let outside_face_index = self.faces.add(Face {
            point_index: Some(point_index),
            starting_edge_index: first_outside_face_edge_index,
        });
        let mut outside_face_edge_index = first_outside_face_edge_index;
        self.edges.get_mut_or_fail(first_outside_face_edge_index).face = Some(outside_face_index);

        // TODO check if it's fine to do this on-demand instead of batching
        let mut vertices_to_destroy = Vec::<Option<usize>>::new();
        let mut edges_to_destroy = Vec::<Option<usize>>::new();

        loop {
            let mut previous_vertex_index = self.edges.get_or_fail(outgoing_edge_index).target;
            vertices_to_destroy.push(previous_vertex_index);

            println!("vertex: {:?}", self.vertices.get_or_fail(previous_vertex_index.unwrap()));

            // The target vertex should not be inside.
            debug_assert_ne!(
                plane.vector_location(
                    self.vertices.get_or_fail(previous_vertex_index.unwrap()),
                    Polyhedron::tolerance()
                ),
                PlaneLocation::Inside
            );

            // Set the next edge and target vertex to potentially cut.
            let mut current_edge_index = self.edges.get_or_fail(outgoing_edge_index).next;
            let mut current_vertex_index = self.target_index(current_edge_index);

            let mut previous_location = plane.vector_location(
                previous_vertex_index
                    .and_then(|i| self.vertices.get(i))
                    .unwrap(),
                Polyhedron::tolerance(),
            );

            // TODO: Is there a better way to check this?
            // TODO: Is there any reason to set this at all?
            // If the starting vertex was outside the plane, we definitely need to cut.
            let mut need_to_cut = previous_location == PlaneLocation::Outside;

            let mut current_location = plane.vector_location(
                current_vertex_index
                    .and_then(|i| self.vertices.get(i))
                    .unwrap(),
                Polyhedron::tolerance(),
            );

            // Keep going until we are back inside.
            while current_location != PlaneLocation::Inside {
                need_to_cut = true;
                vertices_to_destroy.push(current_vertex_index);
                edges_to_destroy.push(current_edge_index);

                previous_vertex_index = current_vertex_index;
                previous_location = current_location;

                current_edge_index = self.next_index(current_edge_index);
                current_vertex_index = self.target_index(current_edge_index);

                current_location = plane.vector_location(
                    self.vertices.get_or_fail(current_vertex_index.unwrap()),
                    Polyhedron::tolerance(),
                );
            }

            println!("v to d {:?}", vertices_to_destroy);
            println!("which is {:?}", self.vertices.get_or_fail(6));

            // TODO: This should cause a bug if `previous_intersection` is not set in the loop.
            self.edges.get_mut_or_fail(outgoing_edge_index).target = previous_intersection;

            if need_to_cut {
                // TODO: Refactor this massive condition into its own function.
                let current_intersection_vertex_index =
                    if previous_location == PlaneLocation::Incident {
                        let old_vertex = self.vertices.get_or_fail(previous_vertex_index.unwrap());

                        // TODO: Can this be cloned instead?
                        let vertex: Vector3<Real> = Vector3 {
                            x: old_vertex.x,
                            y: old_vertex.y,
                            z: old_vertex.z,
                        };

                        self.vertices.add(vertex)
                    } else {
                        let vertex = plane.intersection(
                            self.vertices.get_or_fail(previous_vertex_index.unwrap()),
                            self.vertices.get_or_fail(current_vertex_index.unwrap()),
                        );

                        self.vertices.add(vertex)
                    };

                let current_face_index = self.edges.get_or_fail(outgoing_edge_index).face;

                // Make sure that the face's starting edge is one that won't end up getting destroyed.
                self.faces
                    .get_mut_or_fail(current_face_index.unwrap())
                    .starting_edge_index = outgoing_edge_index;

                let bridge_edge_index = self.edges.add(HalfEdge {
                    face: current_face_index,
                    target: Some(current_intersection_vertex_index),
                    flip: Some(outside_face_edge_index),
                    next: current_edge_index,
                });

                self.edges.get_mut_or_fail(outside_face_edge_index).flip = Some(bridge_edge_index);
                self.edges.get_mut_or_fail(outgoing_edge_index).next = Some(bridge_edge_index);

                outside_face_edge_index = self.edges.add(HalfEdge {
                    face: Some(outside_face_index),
                    target: Some(current_intersection_vertex_index),
                    // TODO: This seems wrong
                    flip: None,
                    next: Some(outside_face_edge_index),
                });

                previous_intersection = Some(current_intersection_vertex_index);
            }

            outgoing_edge_index = self
                .edges
                .get_or_fail(current_edge_index.unwrap())
                .flip
                .unwrap();

            // The test condition that should break out of the loop.
            // println!("first outgoing_edge: {:?}", self.edges.get_or_fail(first_outgoing_edge_index));
            // println!("first outgoing_edge.target: {:?}", self.vertices.get_or_fail(self.edges.get_or_fail(first_outgoing_edge_index).target.unwrap()));
            // println!("outgoing_edge: {:?}", self.edges.get_or_fail(outgoing_edge_index));
            // println!("outgoing_edge.target: {:?}", self.vertices.get_or_fail(self.edges.get_or_fail(outgoing_edge_index).target.unwrap()));
            println!("{:?}", self.edges);
            for e in &self.edges {
                if e.target.is_none() {
                    println!("{:?}", e);
                }
            }
            if outgoing_edge_index == first_outgoing_edge_index {
                break;
            }
        }

        while self
            .target_index(Some(outgoing_edge_index))
            .and_then(|i| self.edges.get(i))
            .is_some()
        {
            self.edges.get_mut_or_fail(outgoing_edge_index).target = previous_intersection;

            outgoing_edge_index = self
                .next_index(Some(outgoing_edge_index))
                .and_then(|i| self.edges.get_or_fail(i).flip)
                .unwrap();
        }

        self.clean_up_vertices(vertices_to_destroy);
        self.clean_up_edges(&mut edges_to_destroy);

        true
    }

    /// Remove the vertices of the polyhedron that have been cut off by plane cutting.
    fn clean_up_vertices(&mut self, vertices_to_destroy: Vec<Option<usize>>) {
        for vertex_index in vertices_to_destroy {
            match vertex_index {
                Some(index) => self.vertices.remove(index),
                _ => { /* Do nothing */ }
            }
        }
    }

    /// Remove the edges of the polyhedron that have been cut off by plane cutting.
    fn clean_up_edges(&mut self, edges_to_destroy: &mut Vec<Option<usize>>) {
        if edges_to_destroy.len() == 0 {
            return;
        }

        // Note that the entire region to destroy should be a single connected component.
        //
        // We begin by deleting a perimiter around the region to destroy, which marks the
        // boundary. At the same time, store any edges found in the inside of the region to be
        // deleted to serve as a starting point for deleting everything else. Any one edge would do,
        // but since we aren't sure which ones are actually part of the perimiter, we just store
        // them all.
        let mut connected_components_to_destroy =
            Vec::<Option<usize>>::with_capacity(edges_to_destroy.len());

        for edge_index in edges_to_destroy {
            // Store related edges in the deletion region.
            // TODO can this be done earlier to save time?
            connected_components_to_destroy.push(self.flip_index(*edge_index));

            match edge_index {
                Some(index) => self.edges.remove(*index),
                _ => { /* Do nothing */ }
            }
        }

        // Search over the region to destroy and do the destroying.
        for edge_index in connected_components_to_destroy {
            // TODO We thought that only a single non-perimeter edge would be enough, since but it
            // seems that that was not the case in the original implementation. It might be worth
            // looking into again, but for now, try them all.
            match edge_index {
                Some(index) => self.mark_sweep(index),
                _ => { /* Do nothing */ }
            }
        }
    }

    /// Delete the portion of the polyhedron connected to the given edge.
    ///
    /// Note that this function should only be used on portions of the polyhedron that have already
    /// been cut off, since calling it on an edge connected to the main portion would delete the
    /// entire polyhedron.
    fn mark_sweep(&mut self, edge_index: usize) {
        let starting_edge = match self.edges.get(edge_index) {
            Some(edge) => edge,
            _ => return,
        };

        let flip = starting_edge.flip;
        let next = starting_edge.next;
        let target = starting_edge.target;
        let face = starting_edge.face;

        match flip {
            Some(f) => self.mark_sweep(f),
            _ => { /* Do nothing */ }
        }

        match next {
            Some(n) => self.mark_sweep(n),
            _ => { /* Do nothing */ }
        }

        match target {
            Some(t) => self.vertices.remove(t),
            _ => { /* Do nothing */ }
        }

        match face {
            Some(f) => self.faces.remove(f),
            _ => { /* Do nothing */ }
        }

        self.edges.remove(edge_index);
    }

    /// Get the target vertex index of the edge at the given edge index, if possible.
    fn target_index(&self, edge_index: Option<usize>) -> Option<usize> {
        edge_index
            .and_then(|i| self.edges.get(i))
            .and_then(|e| e.target)
    }

    /// Get the source vertex index of the edge at the given edge index, if possible.
    fn source_index(&self, edge_index: Option<usize>) -> Option<usize> {
        edge_index
            .and_then(|i| self.edges.get(i))
            .and_then(|e| self.target_index(e.flip))
    }

    /// Get the next edge from the current edge index, if possible.
    fn next_index(&self, edge_index: Option<usize>) -> Option<usize> {
        edge_index
            .and_then(|i| self.edges.get(i))
            .and_then(|e| e.next)
    }

    // Get the flip edge from the current edge index, if possible.
    fn flip_index(&self, edge_index: Option<usize>) -> Option<usize> {
        edge_index
            .and_then(|i| self.edges.get(i))
            .and_then(|e| e.flip)
    }

    /// Check if the polyhedron has been built yet.
    pub fn is_built(&self) -> bool {
        return self.root_edge.is_some();
    }

    /// Get the normal vector for a face.
    ///
    /// The algorithm used here is to take the first vertex of the face and set it as the anchor.
    /// Then move along the vertices in the face until the first one is reached again. At each
    /// vertex, get the vector between it and the anchor. Take the cross product of that vector and
    /// the one for the previous vertex. Finally, add all of cross products to get the weighted
    /// normal.
    ///
    /// What this is doing is getting the normal vector of each triangle that makes up the surface
    /// (using a triangulation where all triangles share the anchor as a vertex) and adding them
    /// together. This method is generally considered to be quite fast.
    pub fn weighted_normal(&self, face_index: usize) -> Vector3<Real> {
        let mut normal = Vector3::<Real> {
            x: 0.into(),
            y: 0.into(),
            z: 0.into(),
        };

        // The target of the starting edge will serve as the anchor.
        let starting_edge_index = self.faces.get_or_fail(face_index).starting_edge_index;
        let starting_edge = self.edges.get_or_fail(starting_edge_index);
        let target_vertex_index = starting_edge.target.unwrap();
        let anchor = self.vertices.get_or_fail(target_vertex_index);

        // Set the first "current vertex". Note that it must be TWO vertices after the anchor,
        // since it is crossed with the previous one.
        let mut current_edge_index = starting_edge.next.unwrap();
        let mut current_edge = self.edges.get_or_fail(current_edge_index);
        let mut current_vector = self.vertices.get_or_fail(current_edge.target.unwrap()) - anchor;
        current_edge_index = current_edge.next.unwrap();
        current_edge = self.edges.get_or_fail(current_edge_index);

        while current_edge_index != starting_edge_index {
            // TODO clone - this can probably be made more elegant in general
            let previous_vector = current_vector.clone();
            current_vector = self.vertices.get_or_fail(current_edge.target.unwrap()) - anchor;
            normal = &normal + &Vector3::cross(&previous_vector, &current_vector);

            current_edge_index = current_edge.next.unwrap();
            current_edge = self.edges.get_or_fail(current_edge_index);
        }

        normal
    }

    /// Compute the face data for each face and save it. Does nothing if the face data already has
    /// been computed.
    fn compute_face_data(&mut self) {
        if !self.face_data.is_empty() {
            return;
        }

        for i in 0..self.faces.len() {
            if self.faces.has(i) {
                self.face_data.push(FaceData {
                    face_index: i,
                    weighted_normal: self.weighted_normal(i),
                })
            }
        }
    }

    /// Compute the volume of the Polyhedron.
    ///
    /// The algorithm is a method derived from the divergence theorem, which is as follows:
    /// 1. For each face, find an arbitrary point on the face and take its dot product with the
    ///    unit normal vector leaving the face.
    /// 2. Multiply the dot product by the area of the face.
    /// 3. Sum the results for each face.
    /// 4. Divide the total by 3.
    ///
    /// Since we compute the weighted normal by volume instead of the unit normal, steps 1 and 2
    /// are combined into a single step.
    pub fn compute_volume(&mut self) -> Real {
        self.compute_face_data();

        let mut volume = Real::from(0);

        for fd in &self.face_data {
            // Get the first vertex of the face.
            let face = self.faces.get_or_fail(fd.face_index);
            let starting_edge = self.edges.get_or_fail(face.starting_edge_index);
            let target_vertex = self.vertices.get_or_fail(starting_edge.target.unwrap());

            volume = volume + Vector3::dot(target_vertex, &fd.weighted_normal);
        }

        // TODO Why divide by 6 instead of 3? It's probably because the weighted normals are
        // doubled.
        volume / (6.into())
    }

    // TODO test
    /// Shift all vertices in the polyhedron by a vector.
    pub fn translate(&mut self, shift: Vector3<Real>) {
        for i in 0..self.vertices.len() {
            // TODO use an actual mutable iterator over the pool, if that's possible
            match self.vertices.get_mut(i) {
                Some(v) => v.add(&shift),
                None => (),
            }
        }
    }

    // TODO test
    /// Get the indices of the points that share a face with this Polyhedron.
    pub fn compute_neighbors(&self) -> Vec<usize> {
        let mut neighbors = Vec::new();

        for face in &self.faces {
            // TODO maybe allow empty point index, but that seems like an error
            // ^ actually it indicates the case where "initial" cube faces are still present
            neighbors.push(face.point_index.unwrap());
        }

        neighbors
    }

    // TODO test
    /// Get the points at the vertices of the Polyhedron.
    pub fn compute_vertices(&self) -> Vec<Vector3<Real>> {
        let mut vertices = Vec::new();

        for vertex in &self.vertices {
            vertices.push(vertex.clone());
        }

        vertices
    }

    // TODO test
    /// Get the vertices for a face.
    pub fn compute_face_vertices(&self, face_index: usize) -> Vec<Vector3<Real>> {
        let mut face_vertices = Vec::new();

        let start = self.faces.get_or_fail(face_index).starting_edge_index;
        let mut current = start.clone();

        loop {
            let edge = self.edges.get_or_fail(current);

            let vertex_index = edge.target.unwrap();
            let vertex = self.vertices.get_or_fail(vertex_index);

            face_vertices.push(vertex.clone());

            // Get the next edge until we loop back around.
            current = edge.next.unwrap();
            if current == start {
                break;
            }
        }

        face_vertices
    }

    // TODO figure out if this is right
    pub fn get_face_point(&self, face_index: usize) -> usize {
        self.faces.get_or_fail(face_index).point_index.unwrap()
    }

    // TODO maybe just make it public? This is probably slower.
    pub fn get_face_indices(&self) -> Vec<usize> {
        let mut face_indices = Vec::new();

        for face in &self.faces {
            face_indices.push(face.point_index.unwrap());
        }

        face_indices
    }
}

#[cfg(test)]
mod tests {
    use super::{HalfEdge, Polyhedron};
    use crate::float::Float64;
    use crate::vector3::{Plane, Vector3};

    fn make_plane(n_x: f64, n_y: f64, n_z: f64, p_x: f64, p_y: f64, p_z: f64) -> Plane<Float64> {
        Plane::build_from_non_unit_normal_and_point(
            Vector3::new(Float64(n_x), Float64(n_y), Float64(n_z)),
            Vector3::new(Float64(p_x), Float64(p_y), Float64(p_z)),
        )
    }

    #[test]
    fn new_test() {
        let p = Polyhedron::new(
            Float64(-3.0),
            Float64(40.0),
            Float64(-0.2),
            Float64(0.0),
            Float64(100.0),
            Float64(-0.1),
        );

        assert_eq!(p.edges.len(), 24);
        assert_eq!(p.vertices.len(), 8);
        assert_eq!(p.faces.len(), 6);
        assert_eq!(p.face_data.len(), 0);
        assert_eq!(p.root_edge.unwrap(), 0);

        for vertex in &p.vertices {
            assert!(vertex.x == Float64(0.0) || vertex.x == Float64(-3.0));
            assert!(vertex.y == Float64(100.0) || vertex.y == Float64(40.0));
            assert!(vertex.z == Float64(-0.1) || vertex.z == Float64(-0.2));
        }
    }

    #[test]
    fn reset_test() {
        let mut p = Polyhedron::new(
            Float64(-3.0),
            Float64(40.0),
            Float64(-0.2),
            Float64(0.0),
            Float64(100.0),
            Float64(-0.1),
        );

        p.reset(
            Float64(-1.0),
            Float64(-1.0),
            Float64(-1.0),
            Float64(1.0),
            Float64(1.0),
            Float64(1.0),
        );

        assert_eq!(p.edges.len(), 24);
        assert_eq!(p.vertices.len(), 8);
        assert_eq!(p.faces.len(), 6);
        assert_eq!(p.face_data.len(), 0);
        assert_eq!(p.root_edge.unwrap(), 0);

        for vertex in &p.vertices {
            assert!(vertex.x == Float64(1.0) || vertex.x == Float64(-1.0));
            assert!(vertex.y == Float64(1.0) || vertex.y == Float64(-1.0));
            assert!(vertex.z == Float64(1.0) || vertex.z == Float64(-1.0));
        }
    }

    #[test]
    fn build_cube_test() {
        let (root_edge, vertices, faces, edges) = Polyhedron::build_cube(
            Float64(-1.0),
            Float64(-1.0),
            Float64(-1.0),
            Float64(1.0),
            Float64(1.0),
            Float64(1.0),
        );

        assert_eq!(root_edge.unwrap(), 0);
        assert_eq!(vertices.len(), 8);
        assert_eq!(faces.len(), 6);
        assert_eq!(edges.len(), 24);

        for e in &edges {
            // The flip of the flip should get back to the original edge.
            assert_eq!(
                edges.get_or_fail(edges.get_or_fail(e.flip.unwrap()).flip.unwrap()),
                e
            );
        }

        // TODO perhaps test more carefully
    }

    #[test]
    fn find_outgoing_edge_no_cut() {
        let poly = Polyhedron::new(
            Float64(-1.0),
            Float64(-1.0),
            Float64(-1.0),
            Float64(1.0),
            Float64(1.0),
            Float64(1.0),
        );

        // This plane lies outside the polyhedron.
        let plane = Plane::build_from_non_unit_normal_and_point(
            Vector3::new(Float64(1.0), Float64(1.0), Float64(1.0)),
            Vector3::new(Float64(2.0), Float64(2.0), Float64(2.0)),
        );

        assert_eq!(poly.find_outgoing_edge(&plane), None);
    }

    #[test]
    fn find_outgoing_edge_with_cut() {
        let poly = Polyhedron::new(
            Float64(-1.0),
            Float64(-1.0),
            Float64(-1.0),
            Float64(1.0),
            Float64(1.0),
            Float64(1.0),
        );

        // This plane should cut the polyhedron, only removing the corner (1, 1, 1).
        let plane = Plane::build_from_non_unit_normal_and_point(
            Vector3::new(Float64(1.0), Float64(1.0), Float64(1.0)),
            Vector3::new(Float64(0.5), Float64(0.5), Float64(0.5)),
        );

        let edge = poly
            .edges
            .get_or_fail(poly.find_outgoing_edge(&plane).unwrap());

        let target_vertex = poly.vertices.get_or_fail(edge.target.unwrap());

        let source_vertex = poly
            .vertices
            .get_or_fail(poly.edges.get_or_fail(edge.flip.unwrap()).target.unwrap());

        let corner = Vector3::new(Float64(1.0), Float64(1.0), Float64(1.0));
        let adjacent_1 = Vector3::new(Float64(-1.0), Float64(1.0), Float64(1.0));
        let adjacent_2 = Vector3::new(Float64(1.0), Float64(-1.0), Float64(1.0));
        let adjacent_3 = Vector3::new(Float64(1.0), Float64(1.0), Float64(-1.0));

        // The target of the edge should be (1, 1, 1).
        assert_eq!(target_vertex, &corner);

        // The source of the edge should be an adjacent vertex.
        assert!(
            source_vertex == &adjacent_1
                || source_vertex == &adjacent_2
                || source_vertex == &adjacent_3
        );
    }

    #[test]
    fn cut_with_plane() {
        let mut poly = Polyhedron::new(
            Float64(-1.0),
            Float64(-1.0),
            Float64(-1.0),
            Float64(1.0),
            Float64(1.0),
            Float64(1.0),
        );

        // The actual value of the point index doesn't matter to the polyhedron, but can be used to
        // find the new face for our testing purposes.
        let point_index = 100;
        let point = Vector3::new(Float64(1.0), Float64(1.0), Float64(1.0));
        let plane = Plane::halfway_from_origin_to(point);

        poly.cut_with_plane(point_index, &plane);

        // TODO test more thoroughly
    }
}
