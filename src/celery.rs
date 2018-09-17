//! Implementation of a basic 3-D cell array

use std::cmp::{min, Ordering};
use std::f64;
use std::vec::Vec;

/// The ideal number of particles in a cell. Used to determine the size of a cell.
// TODO: Special const handling?
// TODO: Experiment with different values
const CELL_DENSITY: f64 = 1.25;

/// Contains a combination of distance information and indices to a cell. Used for sorting cells by
/// distance when searching for neighbors.
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
struct Celery<'a, PointType: ToCeleryPoint + Default + 'a> {
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

    /// The order in which to search through cells.
    search_order: Vec<DistanceIndex>,
}

impl<'a, PointType: ToCeleryPoint + Default> Celery<'a, PointType> {
    /// Create a cell array from a slice of points.
    fn new(pts: &'a [PointType]) -> Celery<PointType> {
        let bounds = CeleryBounds::new(pts);
        let cell_info = CeleryCellInfo::new(pts, &bounds);

        let cells = Celery::get_cells(pts, &bounds, &cell_info);
        let sorted_indices = Celery::get_sorted_indices(pts, &cells);
        let delimiters = Celery::get_delimiters(pts, &cells, &sorted_indices, &bounds, &cell_info);

        Celery::<PointType> {
            bounds: bounds,
            points: pts,
            cells: cells,
            delimiters: delimiters,
            sorted_indices: sorted_indices,
            cell_info: cell_info,
            search_order: Vec::default(), // TO CHANGE
        }
    }

    /// Get the x-index of the cell from an x-coordinate.
    fn get_x_cell_index(x: f64, bounds: &CeleryBounds, cell_info: &CeleryCellInfo) -> usize {
        // TODO: Is this first condition necessary?
        if x >= bounds.x_max {
            cell_info.cells_per_dimension - 1
        } else {
            let index = (x - bounds.x_min) * cell_info.x_inverse_cell_size;
            min(index as usize, cell_info.cells_per_dimension - 1)
        }
    }

    /// Get the y-index of the cell from an y-coordinate.
    fn get_y_cell_index(y: f64, bounds: &CeleryBounds, cell_info: &CeleryCellInfo) -> usize {
        // TODO: Is this first condition necessary?
        if y >= bounds.y_max {
            cell_info.cells_per_dimension - 1
        } else {
            let index = (y - bounds.y_min) * cell_info.y_inverse_cell_size;
            min(index as usize, cell_info.cells_per_dimension - 1)
        }
    }

    /// Get the z-index of the cell from an z-coordinate.
    fn get_z_cell_index(z: f64, bounds: &CeleryBounds, cell_info: &CeleryCellInfo) -> usize {
        // TODO: Is this first condition necessary?
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
            index_map[i] = i;
        }

        index_map.sort_unstable_by(|a, b| cells[*a].cmp(&cells[*b]));

        index_map
    }

    /// Get the delimiters for the sorted cells.
    fn get_delimiters(
        pts: &'a [PointType],
        cells: &Vec<usize>,
        sorted_indices: &Vec<usize>,
        bounds: &CeleryBounds,
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

                    return delimiters
                }
            }

            delimiters.push(last_delimiter + offset);
        }

        // Go one past the end so delimiters can be checked for the last cell.
        delimiters.push(num_points);

        return delimiters
    }
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
}
