// celery.rs --- Implementation of a basic 3-D cell array for The Tessellator
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

use std::cmp::{max, min, Ordering};
use std::f64;
use std::vec::Vec;

/// The ideal number of particles in a cell. Used to determine the size of a cell.
// TODO: Special const handling?
// TODO: Experiment with different values
const CELL_DENSITY: f64 = 1.25;

/// Contains a combination of distance information and indices to a cell. Used for sorting cells by
/// distance when searching for neighbors.
///
/// The distance between two cells is always determined by the relative cell positions, so we can
/// pre-compute that information for each possible cell relationship.
#[derive(Debug, Default)]
struct DistanceIndex {
    distance: f64,
    i: i32,
    j: i32,
    k: i32,
}

impl PartialEq for DistanceIndex {
    fn eq(&self, other: &DistanceIndex) -> bool {
        self.distance == other.distance
    }
}

impl PartialOrd for DistanceIndex {
    fn partial_cmp(&self, other: &DistanceIndex) -> Option<Ordering> {
        self.distance.partial_cmp(&other.distance)
    }
}

impl Eq for DistanceIndex {}

impl Ord for DistanceIndex {
    fn cmp(&self, other: &DistanceIndex) -> Ordering {
        // TODO: I completely flaunt the IEEE standard for floats here, but I need a total ordering
        // to do the search. I might want to revisit this logic, or somehow ensure that having a
        // value of NaN is not possible. It seems this can be done with a wrapper type for f64.
        if self.distance.is_nan() {
            Ordering::Less
        } else if other.distance.is_nan() {
            Ordering::Greater
        } else {
            self.distance.partial_cmp(&other.distance).unwrap()
        }
    }
}

/// A simple implementation of a point in 3-D space.
#[derive(Debug, Default)]
pub struct CeleryPoint {
    x: f64,
    y: f64,
    z: f64,
}

/// A trait to convert some point-containing data into a CeleryPoint so that the Celery can
/// perform arithmetic on it.
pub trait ToCeleryPoint {
    fn to_celery_point(&self) -> CeleryPoint;
}

/// The x, y, and z boundaries for a Celery.
#[derive(Debug, Default)]
struct CeleryBounds {
    /// The x-value minimum boundary of the cell array.
    x_min: f64,
    /// The x-value maximum boundary of the cell array.
    x_max: f64,
    /// The y-value minimum boundary of the cell array.
    y_min: f64,
    /// The y-value maximum boundary of the cell array.
    y_max: f64,
    /// The z-value minimum boundary of the cell array.
    z_min: f64,
    /// The z-value maximum boundary of the cell array.
    z_max: f64,
}

impl CeleryBounds {
    /// Get the bounds for a collection of points.
    fn new<PointType: ToCeleryPoint>(pts: &[PointType]) -> CeleryBounds {
        let mut xmin = f64::MAX;
        let mut ymin = f64::MAX;
        let mut zmin = f64::MAX;

        let mut xmax = f64::MIN;
        let mut ymax = f64::MIN;
        let mut zmax = f64::MIN;

        for p in pts {
            let cp = p.to_celery_point();

            if cp.x < xmin {
                xmin = cp.x
            }
            if cp.x > xmax {
                xmax = cp.x
            }

            if cp.y < ymin {
                ymin = cp.y
            }
            if cp.y > ymax {
                ymax = cp.y
            }

            if cp.z < zmin {
                zmin = cp.z
            }
            if cp.z > zmax {
                zmax = cp.z
            }
        }

        // TODO: Think about adding epsilon-distance between min and max values.

        CeleryBounds {
            x_min: xmin,
            x_max: xmax,
            y_min: ymin,
            y_max: ymax,
            z_min: zmin,
            z_max: zmax,
        }
    }
}

/// Cell data for a cell array.
#[derive(Debug, Default)]
struct CeleryCellInfo {
    /// The size of a cell in the x-dimension.
    x_cell_size: f64,
    /// The size of a cell in the y-dimension.
    y_cell_size: f64,
    /// The size of a cell in the z-dimension.
    z_cell_size: f64,

    // TODO check if this performace optimization is still necessary. It probably is.
    /// The inverse size of a cell in the x-dimension.
    x_inverse_cell_size: f64,
    /// The inverse size of a cell in the y-dimension.
    y_inverse_cell_size: f64,
    /// The inverse size of a cell in the z-dimension.
    z_inverse_cell_size: f64,

    /// The number of cells in each dimension.
    cells_per_dimension: usize,
}

impl CeleryCellInfo {
    /// Get the computed cell info for a collection of points.
    fn new<PointType: ToCeleryPoint>(pts: &[PointType], bounds: &CeleryBounds) -> CeleryCellInfo {
        let num_points: f64 = pts.len() as f64;

        // The number of cell "lenghts" that should span each dimension in the cell array. The
        // cast to `usize` rounds toward zero, so add 1 to round up.
        let cells_per_dimension = (num_points / CELL_DENSITY).cbrt() as usize + 1;

        let x_size = (bounds.x_max - bounds.x_min) / (cells_per_dimension as f64);
        let y_size = (bounds.y_max - bounds.y_min) / (cells_per_dimension as f64);
        let z_size = (bounds.z_max - bounds.z_min) / (cells_per_dimension as f64);

        let x_inv_size = (cells_per_dimension as f64) / (bounds.x_max - bounds.x_min);
        let y_inv_size = (cells_per_dimension as f64) / (bounds.y_max - bounds.y_min);
        let z_inv_size = (cells_per_dimension as f64) / (bounds.z_max - bounds.z_min);

        CeleryCellInfo {
            x_cell_size: x_size,
            y_cell_size: y_size,
            z_cell_size: z_size,
            x_inverse_cell_size: x_inv_size,
            y_inverse_cell_size: y_inv_size,
            z_inverse_cell_size: z_inv_size,
            cells_per_dimension: cells_per_dimension,
        }
    }
}

/// A 3-dimensional cell array, designed for quickly finding nearby points in order of distance.
#[derive(Debug)]
pub struct Celery<'a, PointType: ToCeleryPoint + Default + 'a> {
    /// The points stored in the cell array.
    points: &'a [PointType],
    /// The cell index of each point. The cell index corresponds to the point in `points` with the
    /// same index.
    cells: Vec<usize>,
    /// The delimiters for the sorted cells. The delimiter vector is indexed by cell number, and
    /// the value at that index is the index in `sorted_indices` that corresponds to the first
    /// point in that cell. Every following point in `sorted_indices` is also in that cell, until
    /// the next delimiter is reached. In practice, this makes iterating through all of the points
    /// in a cell easy, since you can use
    ///
    /// `for sorted_index in delimiters[cell]..delimiters[cell+1] { ... }`
    ///
    /// The last element of the delimiters vector is the index after the last index in
    /// `sorted_indices`, so this method even works for the last cell. If the last cell or several
    /// of the last cells are empty, they are given that index as well, so the method used above
    /// would loop over an empty list.
    delimiters: Vec<usize>,
    /// Indices into the points slice, sorted in order of the cell the point is in.
    sorted_indices: Vec<usize>,

    /// The bounds for the space contained in the Celery.
    bounds: CeleryBounds,
    /// Information about the structure of individual cells.
    cell_info: CeleryCellInfo,

    /// The order of relative cell positions to search through.
    search_order: Vec<DistanceIndex>,
}

impl<'a, PointType: ToCeleryPoint + Default> Celery<'a, PointType> {
    /// Create a cell array from a slice of points.
    pub fn new(pts: &'a [PointType]) -> Celery<PointType> {
        let bounds = CeleryBounds::new(pts);
        let cell_info = CeleryCellInfo::new(pts, &bounds);

        let cells = Celery::get_cells(pts, &bounds, &cell_info);
        let sorted_indices = Celery::get_sorted_indices(pts, &cells);
        let delimiters = Celery::<PointType>::get_delimiters(&cells, &sorted_indices, &cell_info);
        let search_order = Celery::<PointType>::get_search_order(&cell_info);

        Celery::<PointType> {
            bounds: bounds,
            points: pts,
            cells: cells,
            delimiters: delimiters,
            sorted_indices: sorted_indices,
            cell_info: cell_info,
            search_order: search_order,
        }
    }

    /// Get the x-index of the cell from an x-coordinate.
    fn get_x_cell_index(x: f64, bounds: &CeleryBounds, cell_info: &CeleryCellInfo) -> usize {
        // This is necessary because of the potential for floating point error, placing the
        // coordinate just outside the bounds.
        if x >= bounds.x_max {
            cell_info.cells_per_dimension - 1
        } else {
            let index = (x - bounds.x_min) * cell_info.x_inverse_cell_size;
            min(index as usize, cell_info.cells_per_dimension - 1)
        }
    }

    /// Get the y-index of the cell from an y-coordinate.
    fn get_y_cell_index(y: f64, bounds: &CeleryBounds, cell_info: &CeleryCellInfo) -> usize {
        // This is necessary because of the potential for floating point error, placing the
        // coordinate just outside the bounds.
        if y >= bounds.y_max {
            cell_info.cells_per_dimension - 1
        } else {
            let index = (y - bounds.y_min) * cell_info.y_inverse_cell_size;
            min(index as usize, cell_info.cells_per_dimension - 1)
        }
    }

    /// Get the z-index of the cell from an z-coordinate.
    fn get_z_cell_index(z: f64, bounds: &CeleryBounds, cell_info: &CeleryCellInfo) -> usize {
        // This is necessary because of the potential for floating point error, placing the
        // coordinate just outside the bounds.
        if z >= bounds.z_max {
            cell_info.cells_per_dimension - 1
        } else {
            let index = (z - bounds.z_min) * cell_info.z_inverse_cell_size;
            min(index as usize, cell_info.cells_per_dimension - 1)
        }
    }

    /// Get a cell from its x, y, and z cell indices.
    fn get_cell_from_indices(x: usize, y: usize, z: usize, cell_info: &CeleryCellInfo) -> usize {
        let cpd = cell_info.cells_per_dimension;
        x * cpd * cpd + y * cpd + z
    }

    /// Get the cell a point belongs in.
    fn get_cell(pt: CeleryPoint, bounds: &CeleryBounds, cell_info: &CeleryCellInfo) -> usize {
        let x_index = Celery::<PointType>::get_x_cell_index(pt.x, bounds, cell_info);
        let y_index = Celery::<PointType>::get_y_cell_index(pt.y, bounds, cell_info);
        let z_index = Celery::<PointType>::get_z_cell_index(pt.z, bounds, cell_info);

        Celery::<PointType>::get_cell_from_indices(x_index, y_index, z_index, cell_info)
    }

    /// Compute the cell that each individual point should be stored in.
    fn get_cells(
        pts: &'a [PointType],
        bounds: &CeleryBounds,
        cell_info: &CeleryCellInfo,
    ) -> Vec<usize> {
        let mut cells = Vec::with_capacity(pts.len());

        for pt in pts {
            let cpt = pt.to_celery_point();
            cells.push(Celery::<PointType>::get_cell(cpt, bounds, cell_info));
        }

        cells
    }

    /// Get a list of the index of each point ordered by cell location.
    fn get_sorted_indices(pts: &'a [PointType], cells: &Vec<usize>) -> Vec<usize> {
        let size = pts.len();
        let mut index_map = Vec::with_capacity(size);

        // Initialize the index map.
        for i in 0..size {
            index_map.push(i);
        }

        index_map.sort_unstable_by(|a, b| cells[*a].cmp(&cells[*b]));

        index_map
    }

    /// Get the delimiters for the sorted cells.
    fn get_delimiters(
        cells: &Vec<usize>,
        sorted_indices: &Vec<usize>,
        cell_info: &CeleryCellInfo,
    ) -> Vec<usize> {
        let num_points = cells.len();
        let cpd = cell_info.cells_per_dimension;
        let last_cell = cpd * cpd * cpd;
        let mut delimiters = Vec::with_capacity(last_cell + 1);

        // The first cell is cell 0.
        delimiters.push(0);

        for i in 0..last_cell {
            let mut offset = 0;
            let last_delimiter = match delimiters.last() {
                None => 0,
                Some(j) => *j,
            };

            // Mark the delimiter once we reach a point that is in a different cell.
            while i == cells[sorted_indices[last_delimiter + offset]] {
                offset += 1;

                // Avoid boing out of bounds by stopping if the last point has been reached before
                // the last cell (i.e. the last cell is empty).
                if last_delimiter + offset == num_points {
                    for _ in i..last_cell {
                        delimiters.push(num_points);
                    }

                    return delimiters;
                }
            }

            delimiters.push(last_delimiter + offset);
        }

        // Go one past the end so delimiters can be checked for the last cell.
        delimiters.push(num_points);

        return delimiters;
    }

    /// Create a sorted array of cells to search in order of distance. Since all cells are the same
    /// size, the index only needs to understand the relative positions of the cells.
    fn get_search_order(cell_info: &CeleryCellInfo) -> Vec<DistanceIndex> {
        let sq = |x: f64| x * x;

        // The square of the Euclidean distance. There's no need to take the square root, since the
        // distance is only used for sorting.
        let distance = |i: i32, j: i32, k: i32| {
            sq((i as f64) * cell_info.x_cell_size)
                + sq((j as f64) * cell_info.y_cell_size)
                + sq((k as f64) * cell_info.z_cell_size)
        };

        // The maximum index for a cell in a single dimension.
        let max_index: i32 = cell_info.cells_per_dimension as i32 - 1;

        let cb = |x: i32| x * x * x;
        let mut search_order = Vec::with_capacity(cb(2 * max_index + 1) as usize);

        // TODO: Maybe find a more elegant way to store this information. The point is that the
        // cell is always closest to itself.
        search_order.push(DistanceIndex {
            distance: -1.0,
            i: 0,
            j: 0,
            k: 0,
        });

        // Index the cells that are offset in the x, y, and z coordinates.
        for i in 0..max_index {
            // TODO: Put the i-only stuff here.
            for j in 0..max_index {
                // TODO: Put the j-only stuff here.
                for k in 0..max_index {
                    let dist = distance(i, j, k);

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: i + 1,
                        j: j + 1,
                        k: k + 1,
                    });

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: i + 1,
                        j: j + 1,
                        k: -k - 1,
                    });

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: i + 1,
                        j: -j - 1,
                        k: k + 1,
                    });

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: -i - 1,
                        j: j + 1,
                        k: k + 1,
                    });

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: i + 1,
                        j: -j - 1,
                        k: -k - 1,
                    });

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: -i - 1,
                        j: j + 1,
                        k: -k - 1,
                    });

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: -i - 1,
                        j: -j - 1,
                        k: k + 1,
                    });

                    search_order.push(DistanceIndex {
                        distance: dist,
                        i: -i - 1,
                        j: -j - 1,
                        k: -k - 1,
                    });
                }
            }
        }

        // TODO: This can actually be put in an earlier loop.
        // Index the cells that are offset in the x and y coordinates.
        for i in 0..max_index {
            for j in 0..max_index {
                let dist = distance(i, j, 0);

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: i + 1,
                    j: j + 1,
                    k: 0,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: i + 1,
                    j: -j - 1,
                    k: 0,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: -i - 1,
                    j: j + 1,
                    k: 0,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: -i - 1,
                    j: -j - 1,
                    k: 0,
                });
            }
        }

        // Index the cells that are offset in the x and z coordinates.
        for i in 0..max_index {
            for k in 0..max_index {
                let dist = distance(i, 0, k);

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: i + 1,
                    j: 0,
                    k: k + 1,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: i + 1,
                    j: 0,
                    k: -k - 1,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: -i - 1,
                    j: 0,
                    k: k + 1,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: -i - 1,
                    j: 0,
                    k: -k - 1,
                });
            }
        }

        // Index the cells that are offset in the y and z coordinates.
        for j in 0..max_index {
            for k in 0..max_index {
                let dist = distance(0, j, k);

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: 0,
                    j: j + 1,
                    k: k + 1,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: 0,
                    j: j + 1,
                    k: -k - 1,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: 0,
                    j: -j - 1,
                    k: k + 1,
                });

                search_order.push(DistanceIndex {
                    distance: dist,
                    i: 0,
                    j: -j - 1,
                    k: -k - 1,
                });
            }
        }

        // TODO: This can actually be put in an earlier loop.
        // Index the cells that are offset only by the x-coordinate.
        for i in 0..max_index {
            let dist = distance(i, 0, 0);

            search_order.push(DistanceIndex {
                distance: dist,
                i: i + 1,
                j: 0,
                k: 0,
            });

            search_order.push(DistanceIndex {
                distance: dist,
                i: -i - 1,
                j: 0,
                k: 0,
            });
        }

        // Index the cells that are offset only by the y-coordinate.
        for j in 0..max_index {
            let dist = distance(0, j, 0);

            search_order.push(DistanceIndex {
                distance: dist,
                i: 0,
                j: j + 1,
                k: 0,
            });

            search_order.push(DistanceIndex {
                distance: dist,
                i: 0,
                j: -j - 1,
                k: 0,
            });
        }

        // Index the cells that are offset only by the z-coordinate.
        for k in 0..max_index {
            let dist = distance(0, 0, k);

            search_order.push(DistanceIndex {
                distance: dist,
                i: 0,
                j: 0,
                k: k + 1,
            });

            search_order.push(DistanceIndex {
                distance: dist,
                i: 0,
                j: 0,
                k: -k - 1,
            });
        }

        // Now, sort the search order list by closest distance.
        search_order.sort_unstable();

        search_order
    }

    /// Get the smaller of two floating point numbers. When in doubt, return the second one.
    fn min_float(a: f64, b: f64) -> f64 {
        if a < b {
            a
        } else {
            b
        }
    }

    /// Get the larger of two floating point numbers. When in doubt, return the second one.
    fn max_float(a: f64, b: f64) -> f64 {
        if a > b {
            a
        } else {
            b
        }
    }

    /// Check whether any part of a cell is within the specified radius of a given point.
    // TODO: This could probably be refactored.
    fn check_cell_in_range(
        &self,
        x: f64,
        y: f64,
        z: f64,
        radius: f64,
        x_index: usize,
        y_index: usize,
        z_index: usize,
    ) -> bool {
        let sq = |x: f64| x * x;

        // The square of the Euclidean distance. There's no need to take the square root, since we
        // can just compare the squared values.
        let distance_sq = |i: i32, j: i32, k: i32| {
            sq(i as f64 * self.cell_info.x_cell_size)
                + sq(j as f64 * self.cell_info.y_cell_size)
                + sq(k as f64 * self.cell_info.z_cell_size)
        };

        let offset = |coord: usize, index: usize| max(0, (coord as i32 - index as i32).abs() - 1);

        distance_sq(
            offset(
                Celery::<PointType>::get_x_cell_index(x, &self.bounds, &self.cell_info),
                x_index,
            ),
            offset(
                Celery::<PointType>::get_y_cell_index(y, &self.bounds, &self.cell_info),
                y_index,
            ),
            offset(
                Celery::<PointType>::get_z_cell_index(z, &self.bounds, &self.cell_info),
                z_index,
            ),
        ) <= radius * radius
    }

    /// Given 3-D coordinates and a radius, finds all cells within the radius of the coordinates.
    /// Returns a vector of indices to those cells.
    pub fn find_cells_in_radius(&self, x: f64, y: f64, z: f64, radius: f64) -> Vec<usize> {
        let x_low = Celery::<PointType>::max_float(x - radius, self.bounds.x_min);
        let y_low = Celery::<PointType>::max_float(y - radius, self.bounds.y_min);
        let z_low = Celery::<PointType>::max_float(z - radius, self.bounds.z_min);

        let x_high = Celery::<PointType>::min_float(x - radius, self.bounds.x_max);
        let y_high = Celery::<PointType>::min_float(y - radius, self.bounds.y_max);
        let z_high = Celery::<PointType>::min_float(z - radius, self.bounds.z_max);

        let x_min_index =
            Celery::<PointType>::get_x_cell_index(x_low, &self.bounds, &self.cell_info);
        let y_min_index =
            Celery::<PointType>::get_y_cell_index(y_low, &self.bounds, &self.cell_info);
        let z_min_index =
            Celery::<PointType>::get_z_cell_index(z_low, &self.bounds, &self.cell_info);

        let x_max_index =
            Celery::<PointType>::get_x_cell_index(x_high, &self.bounds, &self.cell_info);
        let y_max_index =
            Celery::<PointType>::get_y_cell_index(y_high, &self.bounds, &self.cell_info);
        let z_max_index =
            Celery::<PointType>::get_z_cell_index(z_high, &self.bounds, &self.cell_info);

        let mut cell_indices = Vec::new();

        for i in x_min_index..=x_max_index {
            for j in y_min_index..=y_max_index {
                for k in z_min_index..=z_max_index {
                    // Since we search through cells in a grid, some cells, such as corners, might
                    // lie outside the radius.
                    if self.check_cell_in_range(x, y, z, radius, i, j, k) {
                        let cell_index =
                            Celery::<PointType>::get_cell_from_indices(i, j, k, &self.cell_info);
                        cell_indices.push(cell_index);
                    }
                }
            }
        }

        cell_indices
    }
}

/// Used to search outward from a point to incrementally find the neighbors in approximate order of
/// distance.
#[derive(Debug)]
pub struct ExpandingSearch<'b, 'a: 'b, PointType: ToCeleryPoint + Default + 'a> {
    /// The cell array to search over.
    celery: &'b Celery<'a, PointType>,

    /// The last search index that was searched.
    last_search_index: usize,

    /// The x-index of the cell containing the point.
    x_cell_index: usize,
    /// The y-index of the cell containing the point.
    y_cell_index: usize,
    /// The z-index of the cell containing the point.
    z_cell_index: usize,

    /// Keeps track of whether the search is complete.
    done: bool,
}

#[cfg(test)]
mod tests {
    use super::{Celery, CeleryPoint, DistanceIndex, ToCeleryPoint};

    #[derive(Debug, Default)]
    struct TestPoint(f64, f64, f64);

    impl ToCeleryPoint for TestPoint {
        fn to_celery_point(&self) -> CeleryPoint {
            CeleryPoint {
                x: self.0,
                y: self.1,
                z: self.2,
            }
        }
    }

    #[test]
    fn compare_distance_index() {
        let di1 = DistanceIndex {
            distance: 12.0,
            i: 1,
            j: 1,
            k: 1,
        };
        let di2 = DistanceIndex {
            distance: -10.0,
            i: 1,
            j: 1,
            k: 1,
        };
        let di3 = DistanceIndex {
            distance: 12.0,
            i: 2,
            j: 3,
            k: 4,
        };

        assert!(di2 < di1);
        assert!(di2 <= di1);
        assert!(!(di2 == di1));
        assert!(!(di2 >= di1));
        assert!(!(di2 > di1));

        assert!(di1 == di3);
        assert!(di1 <= di3);
        assert!(di1 >= di3);
        assert!(!(di1 > di3));
        assert!(!(di1 < di3));

        assert!(di3 > di2);
        assert!(di3 >= di2);
        assert!(!(di3 == di2));
        assert!(!(di3 <= di2));
        assert!(!(di3 < di2));
    }

    #[test]
    fn create_celery() {
        let pts = vec![
            TestPoint(1.2, 3.4, 8.3),
            TestPoint(4.2, 7.3, 2.7),
            TestPoint(0.3, 1.7, 9.0),
        ];

        Celery::new(&pts);
    }
}
