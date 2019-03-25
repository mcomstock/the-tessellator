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
use crate::polyhedron::Polyhedron;
use crate::vector3::{BoundingBox, Vector3};

/// Contains a vector of particle objects, a spatial data structure, and a polyhedron representing
/// the container.
pub struct Diagram<'a> {
    /// Whether or not the diagram has been readied for Voronoi computations.
    initialized: bool,

    /// The polyhedron that represents the container of the points.
    container_shape: Polyhedron,

    /// The bounding box for the diagram.
    bounding_box: BoundingBox<f64>,

    /// The particles in the diagram.
    particles: Vec<Vector3<f64>>,

    /// A mapping from internal indices to original indices.
    original_indices: Vec<usize>,

    /// A mapping from original indices to internal indices.
    internal_indices: Vec<usize>,

    /// The groups that contain each particle. The group at an index corresponds to the particle at
    /// that index.
    groups: Vec<usize>,

    /// The cell array used to determine the relationships between points.
    cell_array: Celery<'a, Vector3<f64>>,
}

impl<'a> Diagram<'a> {
    /// The number of particle objects in the diagram.
    fn num_particles(&self) -> usize {
        self.particles.len()
    }

    /// Add a particle to the diagram with a specified group. This should only be done before the
    /// diagram is initialized.
    fn add_particle_with_group(&mut self, particle: Vector3<f64>, group: usize) {
        self.bounding_box.adjust_to_contain(&particle);
        self.particles.push(particle);
        self.groups.push(group);
    }
}
