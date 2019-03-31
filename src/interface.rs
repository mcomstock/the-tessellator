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
use crate::float::Float;
use crate::polyhedron::Polyhedron;
use crate::vector3::{BoundingBox, Vector3};

/// Contains a vector of particle objects, a spatial data structure, and a polyhedron representing
/// the container.
#[derive(Debug)]
pub struct Diagram<'a, Real: Float> {
    /// Whether or not the diagram has been readied for Voronoi computations.
    initialized: bool,

    /// The polyhedron that represents the container of the points.
    container_shape: Polyhedron<Real>,

    /// The bounding box for the diagram.
    bounding_box: BoundingBox<Real>,

    /// The particles in the diagram.
    particles: Vec<Vector3<Real>>,

    /// A mapping from internal indices to original indices.
    original_indices: Vec<usize>,

    /// A mapping from original indices to internal indices.
    internal_indices: Vec<usize>,

    /// The groups that contain each particle. The group at an index corresponds to the particle at
    /// that index.
    groups: Vec<usize>,

    /// The cell array used to determine the relationships between points.
    cell_array: Celery<'a, Real, Vector3<Real>>,
}

impl<'a, Real: Float> Diagram<'a, Real> {
    /// The number of particle objects in the diagram.
    fn num_particles(&self) -> usize {
        self.particles.len()
    }

    /// Add a particle to the diagram with a specified group. This should only be done before the
    /// diagram is initialized.
    fn add_particle_with_group(&mut self, particle: Vector3<Real>, group: usize) {
        self.bounding_box.adjust_to_contain(&particle);
        self.particles.push(particle);
        self.groups.push(group);
    }

    /// Initialize the diagram. This should only be called after all the particles have been added.
    fn initialize(&mut self) {
        // Adding padding is probably not necessary, so it will not be done by default. A bounding
        // box with padding can be supplied if desired.

        debug_assert!(!self.initialized);

        self.sort_particles();

        self.initialized = true;
    }

    /// Sort the particles in Morton order. Morton order (or Z-order) is a method of sorting
    /// multi-dimensional data that preserves the locality of data points.
    fn sort_particles(&mut self) {
    }
}























































































































































































































































































































































































































































































































































































































