use crate::r2::point::Point as R2Point;
use crate::r2::rect::Rect as R2Rect;
use crate::s2::cellid::CellID;
use crate::s2::edge_clipping::interpolate_f64;
use crate::s2::edge_crosser::EdgeCrosser;
use crate::s2::edge_crossings::Crossing;
use crate::s2::point::Point;
use crate::s2::shape::Shape;
use crate::s2::shape_index::{CellRelation, ShapeIndex, ShapeIndexCell, ShapeIndexIterator};
use crate::shape::ShapeType;
use std::collections::{HashMap, HashSet};

/// CrossingType specifies the types of edge crossings to be reported.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum CrossingType {
    /// CrossingTypeInterior specifies that only edges whose interiors
    /// cross should be reported.
    Interior,

    /// CrossingTypeAll specifies that all pairs of crossing edges should
    /// be reported, even if the crossing is at a shared vertex.
    All,
}

/// EdgeMap stores a sorted set of edge ids for each shape.
pub type EdgeMap = HashMap<ShapeType, Vec<i32>>;

/// CrossingEdgeQuery is used to find the Edge IDs of Shapes that are crossed by
/// a given edge(s).
///
/// Note that if you need to query many edges, it is more efficient to declare
/// a single CrossingEdgeQuery instance and reuse it.
///
/// If you want to find *all* the pairs of crossing edges, it is more efficient to
/// use the not yet implemented VisitCrossings in shapeutil.
pub struct CrossingEdgeQuery<'a> {
    index: &'a ShapeIndex,

    // Temporary values used while processing a query.
    a: R2Point,
    b: R2Point,
    iter: ShapeIndexIterator<'a>,

    // Candidate cells generated when finding crossings.
    cells: Vec<&'a ShapeIndexCell>,
}

impl<'a> CrossingEdgeQuery<'a> {
    /// Creates a new CrossingEdgeQuery for the given index.
    pub fn new(index: &'a ShapeIndex) -> Self {
        CrossingEdgeQuery {
            index,
            a: R2Point::new(0.0, 0.0),
            b: R2Point::new(0.0, 0.0),
            iter: index.iterator(),
            cells: Vec::new(),
        }
    }

    /// Returns the set of edges of the shape S that intersect the given edge AB.
    /// If the CrossingType is Interior, then only intersections at a point interior to both
    /// edges are reported, while if it is CrossingTypeAll then edges that share a vertex
    /// are also reported.
    pub fn crossings(
        &mut self,
        a: Point,
        b: Point,
        shape: &ShapeType,
        cross_type: CrossingType,
    ) -> Vec<i32> {
        let mut edges = self.candidates(a, b, shape);
        if edges.is_empty() {
            return edges;
        }

        let mut crosser = EdgeCrosser::new(&a, &b);
        let mut out = 0;
        let n = edges.len();

        for i in 0..n {
            let edge = shape.edge(edges[i] as i64);
            let sign = crosser.crossing_sign(&edge.v0, &edge.v1);

            if (cross_type == CrossingType::All
                && (sign == Crossing::Maybe || sign == Crossing::Cross))
                || (cross_type != CrossingType::All && sign == Crossing::Cross)
            {
                edges[out] = edges[i];
                out += 1;
            }
        }

        if out < n {
            edges.truncate(out);
        }

        edges
    }

    /// Returns the set of all edges in the index that intersect the given
    /// edge AB. If crossType is CrossingTypeInterior, then only intersections at a
    /// point interior to both edges are reported, while if it is CrossingTypeAll
    /// then edges that share a vertex are also reported.
    ///
    /// The edges are returned as a mapping from shape to the edges of that shape
    /// that intersect AB. Every returned shape has at least one crossing edge.
    pub fn crossings_edge_map(&mut self, a: Point, b: Point, cross_type: CrossingType) -> EdgeMap {
        let mut edge_map = self.candidates_edge_map(a, b);
        if edge_map.is_empty() {
            return edge_map;
        }

        let mut crosser = EdgeCrosser::new(&a, &b);

        // Filter out shapes that don't have any crossing edges
        edge_map.retain(|shape, edges| {
            let mut out = 0;
            let n = edges.len();

            for i in 0..n {
                let edge = shape.edge(edges[i] as i64);
                let sign = crosser.crossing_sign(&edge.v0, &edge.v1);

                if (cross_type == CrossingType::All
                    && (sign == Crossing::Maybe || sign == Crossing::Cross))
                    || (cross_type != CrossingType::All && sign == Crossing::Cross)
                {
                    edges[out] = edges[i];
                    out += 1;
                }
            }

            if out == 0 {
                false // Remove this shape from the map
            } else {
                if out < n {
                    edges.truncate(out);
                }
                true // Keep this shape in the map
            }
        });

        edge_map
    }

    /// Returns a superset of the edges of the given shape that intersect
    /// the edge AB.
    pub fn candidates(&mut self, a: Point, b: Point, shape: &ShapeType) -> Vec<i32> {
        // For small loops it is faster to use brute force. The threshold below was
        // determined using benchmarks.
        const MAX_BRUTE_FORCE_EDGES: usize = 27;
        let max_edges = shape.num_edges();
        if max_edges <= MAX_BRUTE_FORCE_EDGES as i64 {
            return (0..max_edges as i32).collect();
        }

        // Compute the set of index cells intersected by the query edge.
        self.get_cells_for_edge(a, b);
        if self.cells.is_empty() {
            return Vec::new();
        }

        // Look up the shape ID in the index
        let shape_id = self.index.id_for_shape(shape);

        // Gather all the edges that intersect those cells and sort them.
        let mut edges = Vec::new();
        for cell in &self.cells {
            if let Some(clipped) = cell.find_by_shape_id(shape_id) {
                edges.extend(&clipped.edges);
            }
        }

        if self.cells.len() > 1 {
            edges = Self::unique_ints(edges);
        }

        edges
    }

    /// Returns a map from shapes to a superset of edges for that shape
    /// that intersect the edge AB.
    pub fn candidates_edge_map(&mut self, a: Point, b: Point) -> EdgeMap {
        let _shape = self.index.shape(0).as_ref().cloned();

        // If there are only a few edges then it's faster to use brute force. We
        // only bother with this optimization when there is a single shape.
        if self.index.len() == 1 {
            // Typically this method is called many times, so it is worth checking
            // whether the edge map is empty or already consists of a single entry for
            // this shape, and skip clearing edge map in that case.
            let shape = self.index.shape(0).as_ref().cloned();
            let inner_shape = shape.unwrap();

            // Note that we leave the edge map non-empty even if there are no candidates
            // (i.e., there is a single entry with an empty set of edges).
            let mut edge_map: EdgeMap = HashMap::new();
            edge_map.insert(inner_shape.clone(), self.candidates(a, b, &inner_shape));
            return edge_map;
        }

        let mut edge_map: EdgeMap = HashMap::new();

        // Compute the set of index cells intersected by the query edge.
        self.get_cells_for_edge(a, b);
        if self.cells.is_empty() {
            return edge_map;
        }

        // Gather all the edges that intersect those cells.
        for cell in &self.cells {
            for clipped in &cell.shapes {
                if let Some(s) = self.index.shape(clipped.shape_id) {
                    let edges = edge_map.entry(s).or_insert_with(Vec::new);
                    edges.extend(&clipped.edges);
                }
            }
        }

        // If the edge is long, it may have intersected many cells. We need to deduplicate
        // the edges for each shape.
        if self.cells.len() > 1 {
            for (_, edges) in edge_map.iter_mut() {
                *edges = Self::unique_ints(std::mem::take(edges));
            }
        }

        edge_map
    }

    /// Returns a sorted list of unique integers from the input slice.
    fn unique_ints(input: Vec<i32>) -> Vec<i32> {
        let mut set = HashSet::with_capacity(input.len());
        let mut result = Vec::with_capacity(input.len());

        for &i in &input {
            if set.insert(i) {
                result.push(i);
            }
        }

        result.sort();
        result
    }

    /// Populates the cells field with the set of index cells intersected by an edge AB.
    fn get_cells_for_edge(&mut self, a: Point, b: Point) {
        self.cells.clear();

        let segments = crate::s2::edge_clipping::face_segments(a, b);
        for segment in segments {
            self.a = segment.a;
            self.b = segment.b;

            // Optimization: rather than always starting the recursive subdivision at
            // the top level face cell, instead we start at the smallest S2CellId that
            // contains the edge (the edge root cell). This typically lets us skip
            // quite a few levels of recursion since most edges are short.
            let edge_bound = R2Rect::from_points(&[self.a, self.b]);
            let pcell = crate::s2::padded_cell::PaddedCell::from_cell_id(
                CellID::from_face(segment.face as u64),
                0.0,
            );
            let edge_root = pcell.shrink_to_fit(&edge_bound);

            // Now we need to determine how the edge root cell is related to the cells
            // in the spatial index (cellMap). There are three cases:
            //
            //  1. edgeRoot is an index cell or is contained within an index cell.
            //     In this case we only need to look at the contents of that cell.
            //  2. edgeRoot is subdivided into one or more index cells. In this case
            //     we recursively subdivide to find the cells intersected by AB.
            //  3. edgeRoot does not intersect any index cells. In this case there
            //     is nothing to do.
            let relation = self.iter.locate_cell_id(edge_root);
            if relation == CellRelation::Indexed {
                // edgeRoot is an index cell or is contained by an index cell (case 1).
                if let Some(cell) = self.iter.index_cell() {
                    self.cells.push(cell);
                }
            } else if relation == CellRelation::Subdivided {
                // edgeRoot is subdivided into one or more index cells (case 2). We
                // find the cells intersected by AB using recursive subdivision.
                let pcell = if !edge_root.is_face() {
                    crate::s2::padded_cell::PaddedCell::from_cell_id(edge_root, 0.0)
                } else {
                    pcell
                };

                self.compute_cells_intersected(&pcell, edge_bound);
            }
        }
    }

    /// Computes the index cells intersected by the current edge that are
    /// descendants of pcell and adds them to this query's set of cells.
    fn compute_cells_intersected(
        &mut self,
        pcell: &crate::s2::padded_cell::PaddedCell,
        edge_bound: R2Rect,
    ) {
        self.iter.seek(pcell.id.range_min());
        if self.iter.done() || self.iter.cell_id() > pcell.id.range_max() {
            // The index does not contain pcell or any of its descendants.
            return;
        }

        if self.iter.cell_id() == pcell.id {
            // The index contains this cell exactly.
            if let Some(cell) = self.iter.index_cell() {
                self.cells.push(cell);
            }
            return;
        }

        // Otherwise, split the edge among the four children of pcell.
        let center = pcell.middle().lo();

        if edge_bound.x.hi < center.x {
            // Edge is entirely contained in the two left children.
            self.clip_v_axis(edge_bound, center.y, 0, pcell);
            return;
        } else if edge_bound.x.lo >= center.x {
            // Edge is entirely contained in the two right children.
            self.clip_v_axis(edge_bound, center.y, 1, pcell);
            return;
        }

        let child_bounds = self.split_u_bound(edge_bound.clone(), center.x);
        if edge_bound.y.hi < center.y {
            // Edge is entirely contained in the two lower children.
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, 0, 0),
                child_bounds[0].clone(),
            );
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, 1, 0),
                child_bounds[1].clone(),
            );
        } else if edge_bound.y.lo >= center.y {
            // Edge is entirely contained in the two upper children.
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, 0, 1),
                child_bounds[0].clone(),
            );
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, 1, 1),
                child_bounds[1].clone(),
            );
        } else {
            // The edge bound spans all four children. The edge itself intersects
            // at most three children (since no padding is being used).
            self.clip_v_axis(child_bounds[0].clone(), center.y, 0, pcell);
            self.clip_v_axis(child_bounds[1].clone(), center.y, 1, pcell);
        }
    }

    /// Computes the intersected cells recursively for a given padded cell.
    /// Given either the left (i=0) or right (i=1) side of a padded cell pcell,
    /// determine whether the current edge intersects the lower child, upper child,
    /// or both children, and call compute_cells_intersected recursively on those children.
    fn clip_v_axis(
        &mut self,
        edge_bound: R2Rect,
        center: f64,
        i: usize,
        pcell: &crate::s2::padded_cell::PaddedCell,
    ) {
        if edge_bound.y.hi < center {
            // Edge is entirely contained in the lower child.
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, i as u8, 0),
                edge_bound,
            );
        } else if edge_bound.y.lo >= center {
            // Edge is entirely contained in the upper child.
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, i as u8, 1),
                edge_bound,
            );
        } else {
            // The edge intersects both children.
            let child_bounds = self.split_v_bound(edge_bound, center);
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, i as u8, 0),
                child_bounds[0].clone(),
            );
            self.compute_cells_intersected(
                &crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, i as u8, 1),
                child_bounds[1].clone(),
            );
        }
    }

    /// Returns the bound for two children as a result of splitting the
    /// current edge at the given value U.
    fn split_u_bound(&self, edge_bound: R2Rect, u: f64) -> Vec<R2Rect> {
        let v = edge_bound
            .y
            .clamp_point(interpolate_f64(u, self.a.x, self.b.x, self.a.y, self.b.y));

        // diag indicates which diagonal of the bounding box is spanned by AB:
        // it is 0 if AB has positive slope, and 1 if AB has negative slope.
        let diag = if (self.a.x > self.b.x) != (self.a.y > self.b.y) {
            1
        } else {
            0
        };

        Self::split_bound(edge_bound, 0, diag, u, v)
    }

    /// Returns the bound for two children as a result of splitting the
    /// current edge into two child edges at the given value V.
    fn split_v_bound(&self, edge_bound: R2Rect, v: f64) -> Vec<R2Rect> {
        let u = edge_bound
            .x
            .clamp_point(interpolate_f64(v, self.a.y, self.b.y, self.a.x, self.b.x));

        let diag = if (self.a.x > self.b.x) != (self.a.y > self.b.y) {
            1
        } else {
            0
        };

        Self::split_bound(edge_bound, diag, 0, u, v)
    }

    /// Returns the bounds for the two children as a result of splitting
    /// the current edge into two child edges at the given point (u,v). uEnd and vEnd
    /// indicate which bound endpoints of the first child will be updated.
    fn split_bound(edge_bound: R2Rect, u_end: usize, v_end: usize, u: f64, v: f64) -> Vec<R2Rect> {
        let mut child_bounds = vec![edge_bound.clone(), edge_bound.clone()];

        if u_end == 1 {
            child_bounds[0].x.lo = u;
            child_bounds[1].x.hi = u;
        } else {
            child_bounds[0].x.hi = u;
            child_bounds[1].x.lo = u;
        }

        if v_end == 1 {
            child_bounds[0].y.lo = v;
            child_bounds[1].y.hi = v;
        } else {
            child_bounds[0].y.hi = v;
            child_bounds[1].y.lo = v;
        }

        child_bounds
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::s1::angle::Angle;
    use crate::s2::edge_crosser::next_after;
    use crate::s2::random::one_in;
    use crate::shape::Edge;
    use crate::shape::{Chain, ChainPosition, ReferencePoint};
    use rand::{thread_rng, Rng};
    use std::f64;

    // Struct for testing edge queries
    struct EdgeVectorShape {
        edges: Vec<Edge>,
    }

    impl EdgeVectorShape {
        fn new() -> Self {
            EdgeVectorShape { edges: Vec::new() }
        }

        fn add(&mut self, a: Point, b: Point) {
            self.edges.push(Edge { v0: a, v1: b });
        }
    }

    impl Shape for EdgeVectorShape {
        fn num_edges(&self) -> i64 {
            self.edges.len() as i64
        }

        fn edge(&self, e: i64) -> Edge {
            self.edges[e as usize].clone()
        }

        fn reference_point(&self) -> crate::shape::ReferencePoint {
            ReferencePoint::origin(false)
        }

        fn num_chains(&self) -> i64 {
            self.edges.len() as i64
        }

        fn chain(&self, chain_id: i64) -> Chain {
            Chain {
                start: chain_id,
                length: 1,
            }
        }

        fn chain_edge(&self, chain_id: i64, _offset: i64) -> Edge {
            self.edges[chain_id as usize].clone()
        }

        fn chain_position(&self, edge_id: i64) -> ChainPosition {
            ChainPosition {
                chain_id: edge_id,
                offset: 0,
            }
        }

        fn dimension(&self) -> i64 {
            1 // Dimension of polylines
        }
    }

    // Generate sub-edges of some given edge (a,b).
    // The length of the sub-edges is distributed exponentially over a large range,
    // and the endpoints may be slightly perturbed to one side of (a,b) or the other.
    fn generate_perturbed_sub_edges(a: Point, b: Point, count: usize) -> Vec<Edge> {
        let mut edges = Vec::with_capacity(count);

        let a = a.normalize();
        let b = b.normalize();

        let length0 = a.distance(&b);
        for _ in 0..count {
            let length = length0
                * Angle(f64::powf(
                    1e-15,
                    rand::thread_rng().gen_range(0.0..f64::MAX),
                ));
            let offset = (length0 - length) * Angle(rand::thread_rng().gen_range(0.0..f64::MAX));
            edges.push(Edge {
                v0: perturb_at_distance(offset, a, b),
                v1: perturb_at_distance(offset + length, a, b),
            });
        }

        edges
    }

    fn perturb_at_distance(distance: Angle, a0: Point, b0: Point) -> Point {
        let mut x = crate::s2::edgeutil::interpolate_at_distance(&distance, &a0, &b0);
        let mut r = thread_rng();
        if one_in(&mut r, 2) {
            if one_in(&mut r, 2) {
                x.0.x = next_after(x.0.x, 1.0);
            } else {
                x.0.x = next_after(x.0.x, -1.0);
            }
            if one_in(&mut r, 2) {
                x.0.y = next_after(x.0.y, 1.0);
            } else {
                x.0.y = next_after(x.0.y, -1.0);
            }
            if one_in(&mut r, 2) {
                x.0.z = next_after(x.0.z, 1.0);
            } else {
                x.0.z = next_after(x.0.z, -1.0);
            }

            x = x.normalize();
        }

        x
    }

    // Tests the functionality of the CrossingEdgeQuery
    #[test]
    fn test_unique_ints() {
        let test_cases = vec![
            (vec![], Vec::<i32>::new()),
            (vec![1], vec![1]),
            (vec![3, 2, 1], vec![1, 2, 3]),
            (vec![4, 4, 4], vec![4]),
            (
                vec![
                    1, 2, 3, 4, 2, 3, 5, 4, 6, 1, 2, 3, 4, 5, 7, 8, 1, 3, 1, 2, 3, 9, 3, 2, 1,
                ],
                vec![1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
        ];

        for (input, expected) in test_cases {
            let result = CrossingEdgeQuery::unique_ints(input);
            assert_eq!(result, expected);
        }
    }

    //TODO: Additional tests can be added here - See the remaining tests in `crossing_edge_query_test.go`
}
