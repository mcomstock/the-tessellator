// interface.rs --- The interface classes for Voronoi tessellation using the-tessellator.
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

use crate::celery::{Celery, ExpandingSearch, ToCeleryPoint};
use crate::float::{Float, Particle};
use crate::polyhedron::Polyhedron;
use crate::vector3::{BoundingBox, Plane, Vector3};

/// Contains a vector of particle objects, a spatial data structure, and a polyhedron representing
/// the container.
#[derive(Debug, Default)]
pub struct Diagram<Real: Float, PointType: Particle<Real>> {
    /// Whether or not the diagram has been readied for Voronoi computations.
    initialized: bool,

    /// The polyhedron that represents the container of the points.
    container_shape: Polyhedron<Real>,

    /// The bounding box for the diagram.
    bounding_box: BoundingBox<Real>,

    /// A mapping from internal indices to original indices.
    original_indices: Vec<usize>,

    /// A mapping from original indices to internal indices.
    internal_indices: Vec<usize>,

    /// The groups that contain each particle. The group at an index corresponds to the particle at
    /// that index.
    groups: Vec<usize>,

    /// The cell array used to determine the relationships between points.
    cell_array: Celery<Real, PointType>,
}

impl<Real: Float, PointType: Particle<Real>> Diagram<Real, PointType> {
    /// Add a particle to the diagram with a specified group. This should only be done before the
    /// diagram is initialized.
    fn add_particle_with_group(&mut self, particle: PointType, group: usize) {
        self.bounding_box
            .adjust_to_contain(particle.get_x(), particle.get_y(), particle.get_z());
        self.cell_array.points.push(particle);
        self.groups.push(group);
    }

    /// Initialize the diagram. This should only be called after all the particles have been added.
    fn initialize(&mut self) {
        // Adding padding is probably not necessary, so it will not be done by default. A bounding
        // box with padding can be supplied if desired.

        // Don't allow initializing multiple times.
        debug_assert!(!self.initialized);

        self.sort_particles();

        // If a container polyhedron is not given, we need to create one.
        if self.container_shape.is_built() {
            self.container_shape.reset(
                self.bounding_box.low.get_x(),
                self.bounding_box.low.get_y(),
                self.bounding_box.low.get_z(),
                self.bounding_box.high.get_x(),
                self.bounding_box.high.get_y(),
                self.bounding_box.high.get_z(),
            );
        }

        self.cell_array.reset();

        self.initialized = true;
    }

    /// Sort the particles in Morton order. Morton order (or Z-order) is a method of sorting
    /// multi-dimensional data that preserves the locality of data points.
    fn sort_particles(&mut self) {
        let len = self.cell_array.points.len();

        // Use a dummy value for the resize, since it will get overwritten.
        self.internal_indices.resize(len, 0);
        self.original_indices.resize(len, 0);

        // Initialize the new index arrays.
        for i in 0..len {
            self.internal_indices[i] = i;
            // Map the original particle indices to the Morton index.
            self.original_indices[i] = self.morton_index(i);
        }

        // Sort the internally stored indices.
        let original = &self.original_indices;
        self.internal_indices
            .sort_unstable_by(|a, b| original[*a].cmp(&original[*b]));

        // TODO consider using less of a hack
        let completed = usize::max_value();

        for cycle_start in 0..len {
            // TODO is this actually necessary?
            if self.internal_indices[cycle_start] == completed {
                continue;
            }

            // TODO get rid of clone if possible
            let start_particle = self.cell_array.points[cycle_start].clone();
            let start_group = self.groups[cycle_start];

            let mut current = cycle_start;
            let mut next = self.internal_indices[current];

            // TODO pretty sure the line below isn't necessary because it happens in the loop.
            // self.original_indices[current] = self.internal_indices[current];

            while next != cycle_start {
                self.original_indices[current] = self.internal_indices[current];
                self.internal_indices[current] = completed;

                // TODO look for a better way to deal with this
                self.cell_array.points[current] = self.cell_array.points[next].clone();
                self.groups[current] = self.groups[next];

                current = next;
                next = self.internal_indices[next];
            }

            self.original_indices[current] = self.internal_indices[current];
            self.internal_indices[current] = completed;

            self.cell_array.points[current] = start_particle;
            self.groups[current] = start_group;

            for i in 0..len {
                self.internal_indices[self.original_indices[i]] = i;
            }
        }
    }

    /// Get the Morton index of the particle at the given index.
    fn morton_index(&self, index: usize) -> usize {
        let bounding_box = &self.bounding_box;
        let particle = &self.cell_array.points[index];

        // TODO Why the scaling value? That's what we used, but I don't know why.
        let side_lengths = (&bounding_box.high - &bounding_box.low).scale((1.0 / 1024.0).into());

        // TODO do this in a reasonable way
        let offset = &Vector3::<Real> {
            x: particle.get_x(),
            y: particle.get_y(),
            z: particle.get_z(),
        } - &bounding_box.low;

        let x_index: usize = (offset.x / side_lengths.x).into();
        let y_index: usize = (offset.y / side_lengths.y).into();
        let z_index: usize = (offset.z / side_lengths.z).into();

        let mut morton = 0;
        let mut bit = 0;

        // TODO Why 10 as the "number of levels"? Look up Morton index to figure that out.
        for level in 0..10 {
            morton |= ((x_index >> level) & 1) << bit;
            bit += 1;
            morton |= ((y_index >> level) & 1) << bit;
            bit += 1;
            morton |= ((z_index >> level) & 1) << bit;
            bit += 1;
        }

        morton
    }

    /// Get the cell for the point at an index.
    pub fn get_cell_at_index(
        &self,
        index: usize,
        polyhedron: Polyhedron<Real>,
        search_radius: Option<Real>,
        target_group: Option<usize>,
    ) -> Cell<Real, PointType> {
        let point = &self.cell_array.points[self.internal_indices[index]];
        let position = Vector3 {
            x: point.get_x(),
            y: point.get_y(),
            z: point.get_z(),
        };

        Cell {
            diagram: &self,
            index: Some(index),
            position: position,
            polyhedron: polyhedron,
            search_radius: search_radius,
            target_group: target_group,
        }
    }

    /// Get the cell for the given point, assumed not to be in the diagram.
    pub fn get_cell_at_particle(
        &self,
        point: PointType,
        polyhedron: Polyhedron<Real>,
        search_radius: Option<Real>,
        target_group: Option<usize>,
    ) -> Cell<Real, PointType> {
        let position = Vector3 {
            x: point.get_x(),
            y: point.get_y(),
            z: point.get_z(),
        };

        Cell {
            diagram: &self,
            index: None,
            position: position,
            polyhedron: polyhedron,
            search_radius: search_radius,
            target_group: target_group,
        }
    }
}

/// A Voronoi cell.
#[derive(Debug)]
pub struct Cell<'a, Real: Float, PointType: Particle<Real>> {
    /// The diagram this cell belongs to.
    diagram: &'a Diagram<Real, PointType>,
    /// The index of the point in this cell in the diagram. If not present, then this cell is
    /// centered around a point not in the diagram.
    index: Option<usize>,
    /// The position of the point in the cell.
    position: Vector3<Real>,
    /// The Polyhedron representing the Voronoi cell.
    polyhedron: Polyhedron<Real>,
    /// The search radius for neighbors of the particle in the cell. If this is not provided, an
    /// expanding search will be used instead.
    search_radius: Option<Real>,
    /// The target group of the cell. If specified, the cell should only be cut by points from that
    /// group. If not specified, the cell will be cut by all points.
    target_group: Option<usize>,
}

impl<'a, Real: Float, PointType: Particle<Real>> Cell<'a, Real, PointType> {
    /// Compute the Voronoi cell for a particle using an expanding search strategy.
    pub fn compute_voronoi_cell(&mut self) {
        let polyhedron = &mut self.polyhedron;
        debug_assert!(!polyhedron.is_built());

        // Translate the polyhedron so all the center is the origin. This makes finding the cutting
        // planes easier. It might also have the added benefit of making the floating point
        // arithmetic more accurate.
        //
        // Note that it does not actually affect the points stored in the cell array.
        polyhedron.translate(-&self.position);

        let mut expanding_search = ExpandingSearch::new(
            &self.diagram.cell_array,
            self.position.x,
            self.position.y,
            self.position.z,
        );

        let search_points = match self.search_radius {
            Some(radius) => expanding_search.expand_all_in_radius(radius),
            None => expanding_search.expand_all_no_radius(),
        };

        match (self.target_group, self.index) {
            (Some(target), Some(index)) => {
                for search_point_index in search_points {
                    if search_point_index != index
                        && self.diagram.groups[search_point_index] == target
                    {
                        self.cut_with_point(search_point_index);
                    }
                }
            }

            (Some(target), None) => {
                for search_point_index in search_points {
                    if self.diagram.groups[search_point_index] == target {
                        self.cut_with_point(search_point_index);
                    }
                }
            }

            (None, Some(index)) => {
                for search_point_index in search_points {
                    if search_point_index != index {
                        self.cut_with_point(search_point_index);
                    }
                }
            }

            (None, None) => {
                for search_point_index in search_points {
                    self.cut_with_point(search_point_index);
                }
            }
        }
    }

    /// Cut the cell's polyhedron with the point at the specified index.
    fn cut_with_point(&mut self, point_index: usize) {
        let search_point = &self.diagram.cell_array.points[point_index];

        // TODO implement this as a vector method, if possible
        // Translate the point so that its position is correct relative to the
        // center.
        let search_point_as_vector = &Vector3 {
            x: search_point.get_x(),
            y: search_point.get_y(),
            z: search_point.get_z(),
        } - &self.position;

        self.polyhedron.cut_with_plane(
            self.diagram.original_indices[point_index],
            // Note that this works because we have translated the
            // polyhedron.
            &Plane::halfway_from_origin_to(search_point_as_vector),
        );
    }

    /// Get the volume of the Voronoi cell.
    pub fn compute_volume(&mut self) -> Real {
        self.polyhedron.compute_volume()
    }

    /// Get the neighbors (by index) of the point in the Voronoi cell.
    pub fn compute_neighbors(&self) -> Vec<usize> {
        self.polyhedron.compute_neighbors()
    }

    /// Get the indices of all the points within a radius of the Voronoi cell. If a target group is
    /// specified, only get neighbors in that group.
    pub fn compute_neighbor_cloud(&self, radius: Real, target_group: Option<usize>) -> Vec<usize> {
        let mut expanding_search = ExpandingSearch::new(
            &self.diagram.cell_array,
            self.position.x,
            self.position.y,
            self.position.z,
        );

        let mut neighbors = expanding_search.expand_all_in_radius(radius);

        // Remove points from other groups if a target group is specified.
        match target_group {
            Some(target) => neighbors.retain(|n| self.diagram.groups[*n] == target),
            _ => (),
        };

        neighbors
    }

    /// Get the vertices of the Voronoi cell.
    pub fn compute_vertices(&self) -> Vec<Vector3<Real>> {
        self.polyhedron.compute_vertices()
    }

    /// Get the faces of the Voronoi cell.
    pub fn compute_faces(&self) -> Vec<VoronoiFace<Real, PointType>> {
        let mut faces = Vec::new();

        for face_index in self.polyhedron.get_face_indices() {
            faces.push(VoronoiFace {
                face_index: face_index,
                cell: &self,
            });
        }

        faces
    }

    /// Get the original index of the particle at the center of this cell, as it was passed in.
    pub fn original_index(&self) -> Option<usize> {
        self.index.map(|i| self.diagram.original_indices[i])
    }
}

/// A face of a Voronoi cell.
pub struct VoronoiFace<'a, 'b: 'a, Real: Float, PointType: Particle<Real>> {
    /// The index of the face in the Polyhedron for the cell.
    face_index: usize,

    /// The Voronoi cell that contains this face.
    cell: &'a Cell<'b, Real, PointType>,
}

impl<'a, 'b: 'a, Real: Float, PointType: Particle<Real>> VoronoiFace<'a, 'b, Real, PointType> {
    /// Get the vertices of the Voronoi face.
    pub fn compute_vertices(&self) -> Vec<Vector3<Real>> {
        self.cell.polyhedron.compute_face_vertices(self.face_index)
    }

    /// Get the area of the Voronoi face.
    pub fn compute_area(&self) -> Real {
        Real::from(0.5) * self.cell.polyhedron.weighted_normal(self.face_index).mag()
    }

    /// Get the index of the neighbor cell that cut this cell to make this face.
    pub fn compute_neighbor(&self) -> usize {
        // TODO Is this correct? Seems like this would get this cell.
        self.cell.polyhedron.get_face_point(self.face_index)
    }
}
