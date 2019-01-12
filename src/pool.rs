// pool.rs --- Implementation of the Pool class, which functions as a pool for storing indexed
// data.
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

use std::ops::Index;

/// A single chunk of data to be stored in a pool. It is either active and contains an object, or
/// it is inactive and contains an index pointing to the next inactive chunk (if there is one).
/// As a result, the inactive chunks form a linked list, while the active chunks store data.
#[derive(Debug)]
pub enum PoolChunk<T> {
    Value(T),
    NextIndex(usize),
    End,
}

/// A pool in memory for storing indexed data. It is basically a vector of chunks with an
/// iterator implementation.
#[derive(Debug, Default)]
pub struct Pool<T> {
    /// The vector that stores the actual data.
    data: Vec<PoolChunk<T>>,

    /// The index of the first empty chunk to replace, if there is one. Note that this is not
    /// necessarily the empty chunk with the smallest index.
    first: Option<usize>,
}

impl<T> Index<usize> for Pool<T> {
    type Output = PoolChunk<T>;

    fn index(&self, index: usize) -> &PoolChunk<T> {
        &self.data[index]
    }
}

impl<'a, T> IntoIterator for &'a Pool<T> {
    type Item = &'a T;
    type IntoIter = PoolIterator<'a, T>;

    fn into_iter(self) -> PoolIterator<'a, T> {
        PoolIterator {
            pool: &self,
            current_index: 0,
            max_index: self.data.len(),
        }
    }
}

impl<T> Pool<T> {
    /// Create a new pool.
    pub fn new() -> Pool<T> {
        Pool::<T> {
            data: Vec::<PoolChunk<T>>::new(),
            first: Option::None,
        }
    }

    /// Create a new pool with reserved memory capacity.
    pub fn with_capacity(capacity: usize) -> Pool<T> {
        Pool::<T> {
            data: Vec::<PoolChunk<T>>::with_capacity(capacity),
            first: Option::None,
        }
    }

    /// Add a value to the pool, and return the index where the value was added.
    // TODO: The pool in the C++ version would construct items in-place in memory. Unfortunately,
    // I don't think there's a way to do that in Rust, although apparenty it might be added soon.
    // It's worth testing this to see if it's possible to take advantage of LLVM optimizations
    // here.
    pub fn add(&mut self, value: T) -> usize {
        match self.first {
            // If there is an empty spot in the pool, use that.
            Some(i) => {
                // Update the first avaliable index.
                self.first = match self.data[i] {
                    PoolChunk::NextIndex(j) => Option::Some(j),
                    PoolChunk::End => Option::None,
                    PoolChunk::Value(_) => {
                        // There should never be a value here.
                        debug_assert!(false);
                        Option::None
                    }
                };

                self.data[i] = PoolChunk::Value(value);
                return i;
            }

            // If the pool is full, just add to the end.
            None => {
                self.data.push(PoolChunk::Value(value));
                return self.data.len() - 1;
            }
        }
    }

    /// Remove the value at an index from the pool.
    pub fn remove(&mut self, index: usize) {
        let chunk = match self.first {
            Some(i) => PoolChunk::NextIndex(i),
            None => PoolChunk::End,
        };

        self.first = Option::Some(index);
        self.data[index] = chunk;
    }
}

/// An iterator for a pool, which skips over any empty values.
#[derive(Debug)]
pub struct PoolIterator<'a, T: 'a> {
    pool: &'a Pool<T>,
    current_index: usize,
    max_index: usize,
}

impl<'a, T> Iterator for PoolIterator<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<&'a T> {
        while self.current_index < self.max_index {
            let value = &self.pool[self.current_index];

            self.current_index += 1;

            match value {
                PoolChunk::Value(t) => return Some(t),
                _ => (),
            }
        }

        return None;
    }
}

#[cfg(test)]
mod tests {
    use super::{
        Pool,
        PoolChunk::{End, NextIndex, Value},
    };

    struct TestStruct(String);

    #[test]
    fn initial_state() {
        let p = Pool::<TestStruct>::new();

        assert_eq!(p.data.len(), 0);
        assert!(p.first.is_none());
    }

    #[test]
    fn add_one() {
        let mut p = Pool::<TestStruct>::new();
        let ts = TestStruct("hello".to_string());

        let i = p.add(ts);

        assert_eq!(p.data.len(), 1);
        assert!(p.first.is_none());
        assert_eq!(i, 0);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "hello".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn add_one_outlive_initialization() {
        let mut p = Pool::<TestStruct>::new();

        {
            let ts = TestStruct("hello".to_string());
            p.add(ts);
        }

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "hello".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn remove_one_from_one() {
        let mut p = Pool::<TestStruct>::new();
        let ts = TestStruct("hello".to_string());

        p.add(ts);
        p.remove(0);

        assert_eq!(p.data.len(), 1);

        match p.first {
            Some(v) => assert_eq!(v, 0),
            _ => assert!(false),
        };

        match p[0] {
            End => (),
            _ => assert!(false),
        };
    }

    #[test]
    fn add_five() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let i0 = p.add(ts0);
        let i1 = p.add(ts1);
        let i2 = p.add(ts2);
        let i3 = p.add(ts3);
        let i4 = p.add(ts4);

        assert_eq!(p.data.len(), 5);
        assert!(p.first.is_none());

        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
        assert_eq!(i2, 2);
        assert_eq!(i3, 3);
        assert_eq!(i4, 4);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "0".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "1".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            Value(ref v) => assert_eq!(v.0, "2".to_string()),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "3".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            Value(ref v) => assert_eq!(v.0, "4".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn remove_one_from_five() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);

        assert_eq!(p.data.len(), 5);

        match p.first {
            Some(v) => assert_eq!(v, 2),
            _ => assert!(false),
        };

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "0".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "1".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            End => (),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "3".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            Value(ref v) => assert_eq!(v.0, "4".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn remove_two_from_five() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        p.remove(4);

        assert_eq!(p.data.len(), 5);

        match p.first {
            Some(v) => assert_eq!(v, 4),
            _ => assert!(false),
        };

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "0".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "1".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            End => (),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "3".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            NextIndex(i) => assert_eq!(i, 2),
            _ => assert!(false),
        };
    }

    #[test]
    fn remove_two_from_five_reverse() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(4);
        p.remove(2);

        assert_eq!(p.data.len(), 5);

        match p.first {
            Some(v) => assert_eq!(v, 2),
            _ => assert!(false),
        };

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "0".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "1".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            NextIndex(i) => assert_eq!(i, 4),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "3".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            End => (),
            _ => assert!(false),
        };
    }

    #[test]
    fn remove_all_from_five() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(0);
        p.remove(1);
        p.remove(2);
        p.remove(3);
        p.remove(4);

        assert_eq!(p.data.len(), 5);

        match p.first {
            Some(v) => assert_eq!(v, 4),
            _ => assert!(false),
        };

        match p[0] {
            End => (),
            _ => assert!(false),
        };

        match p[1] {
            NextIndex(i) => assert_eq!(i, 0),
            _ => assert!(false),
        };

        match p[2] {
            NextIndex(i) => assert_eq!(i, 1),
            _ => assert!(false),
        };

        match p[3] {
            NextIndex(i) => assert_eq!(i, 2),
            _ => assert!(false),
        };

        match p[4] {
            NextIndex(i) => assert_eq!(i, 3),
            _ => assert!(false),
        };
    }

    #[test]
    fn remove_all_from_five_reverse() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(4);
        p.remove(3);
        p.remove(2);
        p.remove(1);
        p.remove(0);

        assert_eq!(p.data.len(), 5);

        match p.first {
            Some(v) => assert_eq!(v, 0),
            _ => assert!(false),
        };

        match p[0] {
            NextIndex(i) => assert_eq!(i, 1),
            _ => assert!(false),
        };

        match p[1] {
            NextIndex(i) => assert_eq!(i, 2),
            _ => assert!(false),
        };

        match p[2] {
            NextIndex(i) => assert_eq!(i, 3),
            _ => assert!(false),
        };

        match p[3] {
            NextIndex(i) => assert_eq!(i, 4),
            _ => assert!(false),
        };

        match p[4] {
            End => (),
            _ => assert!(false),
        };
    }

    #[test]
    fn replace_one_from_one() {
        let mut p = Pool::<TestStruct>::new();

        let ts = TestStruct("hello".to_string());
        let ts_replace = TestStruct("there".to_string());

        p.add(ts);
        p.remove(0);
        let i0 = p.add(ts_replace);

        assert_eq!(p.data.len(), 1);
        assert!(p.first.is_none());
        assert_eq!(i0, 0);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "there".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn replace_one_from_five() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts2_replace = TestStruct("22".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        let i2 = p.add(ts2_replace);

        assert_eq!(p.data.len(), 5);
        assert!(p.first.is_none());
        assert_eq!(i2, 2);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "0".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "1".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            Value(ref v) => assert_eq!(v.0, "22".to_string()),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "3".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            Value(ref v) => assert_eq!(v.0, "4".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn replace_two_from_five() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts2_replace = TestStruct("22".to_string());
        let ts4_replace = TestStruct("44".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        p.remove(4);

        let i4 = p.add(ts4_replace);
        let i2 = p.add(ts2_replace);

        assert_eq!(p.data.len(), 5);
        assert!(p.first.is_none());

        assert_eq!(i4, 4);
        assert_eq!(i2, 2);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "0".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "1".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            Value(ref v) => assert_eq!(v.0, "22".to_string()),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "3".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            Value(ref v) => assert_eq!(v.0, "44".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn replace_two_from_five_reverse() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts2_replace = TestStruct("22".to_string());
        let ts4_replace = TestStruct("44".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(4);
        p.remove(2);

        let i2 = p.add(ts2_replace);
        let i4 = p.add(ts4_replace);

        assert_eq!(p.data.len(), 5);
        assert!(p.first.is_none());

        assert_eq!(i2, 2);
        assert_eq!(i4, 4);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "0".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "1".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            Value(ref v) => assert_eq!(v.0, "22".to_string()),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "3".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            Value(ref v) => assert_eq!(v.0, "44".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn replace_all_from_five() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts0_replace = TestStruct("00".to_string());
        let ts1_replace = TestStruct("11".to_string());
        let ts2_replace = TestStruct("22".to_string());
        let ts3_replace = TestStruct("33".to_string());
        let ts4_replace = TestStruct("44".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(0);
        p.remove(1);
        p.remove(2);
        p.remove(3);
        p.remove(4);

        let i4 = p.add(ts4_replace);
        let i3 = p.add(ts3_replace);
        let i2 = p.add(ts2_replace);
        let i1 = p.add(ts1_replace);
        let i0 = p.add(ts0_replace);

        assert_eq!(p.data.len(), 5);
        assert!(p.first.is_none());

        assert_eq!(i4, 4);
        assert_eq!(i3, 3);
        assert_eq!(i2, 2);
        assert_eq!(i1, 1);
        assert_eq!(i0, 0);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "00".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "11".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            Value(ref v) => assert_eq!(v.0, "22".to_string()),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "33".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            Value(ref v) => assert_eq!(v.0, "44".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn replace_all_from_five_reverse() {
        let mut p = Pool::<TestStruct>::with_capacity(5);

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts0_replace = TestStruct("00".to_string());
        let ts1_replace = TestStruct("11".to_string());
        let ts2_replace = TestStruct("22".to_string());
        let ts3_replace = TestStruct("33".to_string());
        let ts4_replace = TestStruct("44".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(4);
        p.remove(3);
        p.remove(2);
        p.remove(1);
        p.remove(0);

        let i0 = p.add(ts0_replace);
        let i1 = p.add(ts1_replace);
        let i2 = p.add(ts2_replace);
        let i3 = p.add(ts3_replace);
        let i4 = p.add(ts4_replace);

        assert_eq!(p.data.len(), 5);
        assert!(p.first.is_none());

        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
        assert_eq!(i2, 2);
        assert_eq!(i3, 3);
        assert_eq!(i4, 4);

        match p[0] {
            Value(ref v) => assert_eq!(v.0, "00".to_string()),
            _ => assert!(false),
        };

        match p[1] {
            Value(ref v) => assert_eq!(v.0, "11".to_string()),
            _ => assert!(false),
        };

        match p[2] {
            Value(ref v) => assert_eq!(v.0, "22".to_string()),
            _ => assert!(false),
        };

        match p[3] {
            Value(ref v) => assert_eq!(v.0, "33".to_string()),
            _ => assert!(false),
        };

        match p[4] {
            Value(ref v) => assert_eq!(v.0, "44".to_string()),
            _ => assert!(false),
        };
    }

    #[test]
    fn loop_empty() {
        let p = Pool::<TestStruct>::new();

        for _ in p.into_iter() {
            assert!(false);
        }
    }

    #[test]
    fn loop_one() {
        let mut p = Pool::<TestStruct>::new();
        let ts = TestStruct("hello".to_string());

        p.add(ts);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "hello"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 1);
    }

    #[test]
    fn loop_five() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "0"),
                1 => assert_eq!(v.0, "1"),
                2 => assert_eq!(v.0, "2"),
                3 => assert_eq!(v.0, "3"),
                4 => assert_eq!(v.0, "4"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 5);
    }

    #[test]
    fn loop_one_removed() {
        let mut p = Pool::<TestStruct>::new();
        let ts = TestStruct("hello".to_string());

        p.add(ts);
        p.remove(0);

        for _ in p.into_iter() {
            assert!(false);
        }
    }

    #[test]
    fn loop_five_one_removed() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "0"),
                1 => assert_eq!(v.0, "1"),
                2 => assert_eq!(v.0, "3"),
                3 => assert_eq!(v.0, "4"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 4);
    }

    #[test]
    fn loop_five_two_removed() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        p.remove(4);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "0"),
                1 => assert_eq!(v.0, "1"),
                2 => assert_eq!(v.0, "3"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 3);
    }

    #[test]
    fn loop_five_all_removed() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        p.remove(4);
        p.remove(0);
        p.remove(1);
        p.remove(3);

        for _ in p.into_iter() {
            assert!(false);
        }
    }

    #[test]
    fn loop_five_one_replaced() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts2_replace = TestStruct("22".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        p.add(ts2_replace);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "0"),
                1 => assert_eq!(v.0, "1"),
                2 => assert_eq!(v.0, "22"),
                3 => assert_eq!(v.0, "3"),
                4 => assert_eq!(v.0, "4"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 5);
    }

    #[test]
    fn loop_five_two_replaced() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts2_replace = TestStruct("22".to_string());
        let ts4_replace = TestStruct("44".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        p.remove(4);

        p.add(ts4_replace);
        p.add(ts2_replace);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "0"),
                1 => assert_eq!(v.0, "1"),
                2 => assert_eq!(v.0, "22"),
                3 => assert_eq!(v.0, "3"),
                4 => assert_eq!(v.0, "44"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 5);
    }

    #[test]
    fn loop_five_all_replaced() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts0_replace = TestStruct("00".to_string());
        let ts1_replace = TestStruct("11".to_string());
        let ts2_replace = TestStruct("22".to_string());
        let ts3_replace = TestStruct("33".to_string());
        let ts4_replace = TestStruct("44".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(0);
        p.remove(1);
        p.remove(2);
        p.remove(3);
        p.remove(4);

        p.add(ts4_replace);
        p.add(ts3_replace);
        p.add(ts2_replace);
        p.add(ts1_replace);
        p.add(ts0_replace);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "00"),
                1 => assert_eq!(v.0, "11"),
                2 => assert_eq!(v.0, "22"),
                3 => assert_eq!(v.0, "33"),
                4 => assert_eq!(v.0, "44"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 5);
    }

    #[test]
    fn loop_five_one_removed_one_replaced() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts4_replace = TestStruct("44".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(2);
        p.remove(4);

        p.add(ts4_replace);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "0"),
                1 => assert_eq!(v.0, "1"),
                2 => assert_eq!(v.0, "3"),
                3 => assert_eq!(v.0, "44"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 4);
    }

    #[test]
    fn loop_five_three_removed_two_replaced() {
        let mut p = Pool::<TestStruct>::new();

        let ts0 = TestStruct("0".to_string());
        let ts1 = TestStruct("1".to_string());
        let ts2 = TestStruct("2".to_string());
        let ts3 = TestStruct("3".to_string());
        let ts4 = TestStruct("4".to_string());

        let ts1_replace = TestStruct("11".to_string());
        let ts3_replace = TestStruct("33".to_string());

        p.add(ts0);
        p.add(ts1);
        p.add(ts2);
        p.add(ts3);
        p.add(ts4);

        p.remove(0);
        p.remove(2);
        p.remove(4);
        p.remove(1);
        p.remove(3);

        p.add(ts3_replace);
        p.add(ts1_replace);

        let mut num_iter = 0;
        for (i, v) in p.into_iter().enumerate() {
            match i {
                0 => assert_eq!(v.0, "11"),
                1 => assert_eq!(v.0, "33"),
                _ => assert!(false),
            };

            num_iter = i + 1;
        }

        assert_eq!(num_iter, 2);
    }
}
