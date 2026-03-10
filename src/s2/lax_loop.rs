// Copyright 2023 Google Inc. All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

use crate::s2::point::Point;
use crate::shape::{Chain, ChainPosition, Edge, ReferencePoint, Shape};

use super::rect_bounder::RectBounder;
use super::region::Region;

/// LaxLoop represents a closed loop of edges surrounding an interior
/// region. It is similar to Loop except that this class allows
/// duplicate vertices and edges. Loops may have any number of vertices,
/// including 0, 1, or 2. (A one-vertex loop defines a degenerate edge
/// consisting of a single point.)
///
/// Note that LaxLoop is faster to initialize and more compact than
/// Loop, but does not support the same operations as Loop.
#[derive(Debug, Clone, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LaxLoop {
    pub vertices: Vec<Point>,
}

impl LaxLoop {
    /// Creates a LaxLoop from the given points.
    pub fn from_points(vertices: Vec<Point>) -> Self {
        LaxLoop { vertices }
    }

    /// Creates a LaxLoop from the given Loop, copying its points.
    /// NOTE: This will be implemented once we port Loop to Rust
    // pub fn from_loop(_loop: &Loop) -> Self {
    //     unimplemented!("Requires Loop implementation")
    // }

    /// Returns the vertex at the specified index.
    pub fn vertex(&self, i: usize) -> Point {
        self.vertices[i]
    }
}

impl Region for LaxLoop {
    fn cap_bound(&self) -> super::cap::Cap {
        self.rect_bound().cap_bound()
    }

    fn rect_bound(&self) -> super::rect::Rect {
        let mut bounder = RectBounder::new();
        for v in self.vertices.iter() {
            bounder.add_point(v);
        }
        bounder.get_bound()
    }
}

impl Shape for LaxLoop {
    fn num_edges(&self) -> i64 {
        self.vertices.len() as i64
    }

    fn edge(&self, e: i64) -> Edge {
        let e1 = (e + 1) % self.vertices.len() as i64;
        Edge {
            v0: self.vertices[e as usize],
            v1: self.vertices[e1 as usize],
        }
    }

    fn reference_point(&self) -> ReferencePoint {
        // This is a placeholder that will need a proper implementation
        // In Go, it calls referencePointForShape(l)
        ReferencePoint {
            point: Point::default(),
            contained: false,
        }
    }

    fn num_chains(&self) -> i64 {
        std::cmp::min(1, self.vertices.len() as i64)
    }

    fn chain(&self, _i: i64) -> Chain {
        Chain {
            start: 0,
            length: self.vertices.len() as i64,
        }
    }

    fn chain_edge(&self, _chain_id: i64, offset: i64) -> Edge {
        let j = offset;
        let k = (j + 1) % self.vertices.len() as i64;
        Edge {
            v0: self.vertices[j as usize],
            v1: self.vertices[k as usize],
        }
    }

    fn chain_position(&self, edge_id: i64) -> ChainPosition {
        ChainPosition {
            chain_id: 0,
            offset: edge_id,
        }
    }

    fn dimension(&self) -> i64 {
        2 // Polygons have dimension 2
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lax_loop_from_points() {
        let vertices = vec![
            Point::from_coords(1.0, 0.0, 0.0),
            Point::from_coords(0.0, 1.0, 0.0),
            Point::from_coords(0.0, 0.0, 1.0),
        ];
        let lax_loop = LaxLoop::from_points(vertices.clone());

        assert_eq!(lax_loop.vertices.len(), 3);
        assert_eq!(lax_loop.num_edges(), 3);

        for i in 0..3 {
            assert_eq!(lax_loop.vertex(i), vertices[i]);
        }
    }

    #[test]
    fn test_lax_loop_edge() {
        let vertices = vec![
            Point::from_coords(1.0, 0.0, 0.0),
            Point::from_coords(0.0, 1.0, 0.0),
            Point::from_coords(0.0, 0.0, 1.0),
        ];
        let lax_loop = LaxLoop::from_points(vertices.clone());

        // Test the edges
        for i in 0..3 {
            let edge = lax_loop.edge(i);
            assert_eq!(edge.v0, vertices[i as usize]);
            assert_eq!(edge.v1, vertices[((i + 1) % 3) as usize]);
        }
    }

    #[test]
    fn test_lax_loop_chains() {
        let vertices = vec![
            Point::from_coords(1.0, 0.0, 0.0),
            Point::from_coords(0.0, 1.0, 0.0),
            Point::from_coords(0.0, 0.0, 1.0),
        ];
        let lax_loop = LaxLoop::from_points(vertices);

        assert_eq!(lax_loop.num_chains(), 1);

        let chain = lax_loop.chain(0);
        assert_eq!(chain.start, 0);
        assert_eq!(chain.length, 3);

        // Test chain position
        let pos = lax_loop.chain_position(1);
        assert_eq!(pos.chain_id, 0);
        assert_eq!(pos.offset, 1);
    }

    #[test]
    fn test_empty_lax_loop() {
        let lax_loop = LaxLoop::from_points(vec![]);

        assert_eq!(lax_loop.vertices.len(), 0);
        assert_eq!(lax_loop.num_edges(), 0);
        assert_eq!(lax_loop.num_chains(), 0);
    }
}
