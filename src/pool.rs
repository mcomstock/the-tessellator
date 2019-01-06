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
#[derive(Debug,Default)]
pub struct Pool<T> {
    /// The vector that stores the actual data.
    data: Vec<PoolChunk<T>>,

    /// The index of the first avaliable value, if there is one.
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
            current_index: match self.first {
                Some(val) => val,
                None => 0,
            },
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
    pub fn add(&mut self, value: T) -> usize {
        match self.first {
            // If there is an empty spot in the pool, use that.
            Some(i) => {
                // Update the first avaliable index.
                self.first = match self.data[i] {
                    PoolChunk::NextIndex(j) => Option::Some(j),
                    // Note that there should never be a value in this case.
                    // TODO: Add a debug check for that.
                    _ => Option::None,
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

            match value {
                PoolChunk::Value(t) => return Some(t),
                _ => self.current_index += 1,
            }
        }

        return None;
    }
}
