//! Implementation of a basic 3-D cell array

use std::cmp::Ordering;

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

#[cfg(test)]
mod tests {
    use super::DistanceIndex;

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
}
