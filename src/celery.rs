//! Implementation of a basic 3-D cell array

use std::cmp::Ordering;
use std::vec::Vec;

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

/// A 3-dimensional cell array, designed for quickly finding nearby points in order of distance.
#[derive(Debug, Default)]
struct Celery<PointType: ToCeleryPoint + Default> {
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

    /// The array of points stored in the cell array.
    points: Vec<PointType>,

    /// The cell index of each point. The cell index corresponds to the point in `points` with the
    /// same index.
    cells: Vec<usize>,

    /// The delimiters for each cell.
    delimiters: Vec<usize>,

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
    num_cells: usize,

    /// The order in which to search through cells.
    search_order: Vec<DistanceIndex>,
}

impl<PointType: ToCeleryPoint + Default> Celery<PointType> {
    /// The ideal number of particles in a cell. Used to determine the size of a cell.
    // TODO: Special const handling?
    // TODO: Experiment with different values
    const CELL_DENSITY: f64 = 5.0/4.0;
}

#[cfg(test)]
mod tests {
    use super::{Celery, CeleryPoint, DistanceIndex, ToCeleryPoint};

    #[derive(Debug, Default)]
    struct TestPoint(f64, f64, f64);

    impl ToCeleryPoint for TestPoint {
        fn to_celery_point(&self) -> CeleryPoint {
            CeleryPoint { x: self.0, y: self.1, z: self.2 }
        }
    }

    #[test]
    fn compare_distance_index() {
        let di1 = DistanceIndex { distance: 12.0, i: 1, j: 1, k: 1 };
        let di2 = DistanceIndex { distance: -10.0, i: 1, j: 1, k: 1 };
        let di3 = DistanceIndex { distance: 12.0, i: 2, j: 3, k: 4 };

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
    fn celery_get_cell_density() {
        assert_eq!(Celery::<TestPoint>::CELL_DENSITY, 5.0/4.0);
    }

    #[test]
    fn create_default_celery() {
        Celery::<TestPoint>::default();
    }
}
