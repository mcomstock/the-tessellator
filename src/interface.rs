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

use crate::celery::Celery;
use crate::float::{Float, Particle};
use crate::polyhedron::Polyhedron;
use crate::vector3::{BoundingBox, Vector3};

/// Contains a vector of particle objects, a spatial data structure, and a polyhedron representing
/// the container.
#[derive(Debug)]
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
    /// The number of particle objects in the diagram.
    fn num_particles(&self) -> usize {
        self.cell_array.points.len()
    }

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

        debug_assert!(!self.initialized);

        self.sort_particles();

        // TODO if a container polyhedron is not given, we need to create one.

        self.initialized = true;
    }

    /// Sort the particles in Morton order. Morton order (or Z-order) is a method of sorting
    /// multi-dimensional data that preserves the locality of data points.
    fn sort_particles(&mut self) {
        let len = self.num_particles();

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
}
