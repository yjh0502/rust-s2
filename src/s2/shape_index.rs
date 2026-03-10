use crate::consts::DBL_EPSILON;
use crate::metric::AVG_EDGEMETRIC;
use crate::r1;
use crate::r2::point::Point as R2Point;
use crate::r2::rect::Rect as R2Rect;
use crate::s2::cellid::CellID;
use crate::s2::edge_clipping::interpolate_f64;
use crate::s2::edge_crosser::EdgeCrosser;
use crate::s2::point::Point;
use crate::s2::shape_index_region::{ContainsPointQuery, ShapeIndexRegion, VertexModel};
use crate::s2::stuv::{face, face_uv_to_xyz, valid_face_xyz_to_uv};
use crate::shape::{Edge, Shape, ShapeType};
use crate::shape_index::Status::Fresh;
use std::collections::{BTreeMap, HashMap};
use std::fmt::{Debug, Formatter};
use std::ops::Sub;
use std::ptr::write;
use std::sync::{Arc, RwLock};

// edgeClipErrorUVCoord is the maximum error in a u- or v-coordinate
// compared to the exact result, assuming that the points A and B are in
// the rectangle [-1,1]x[1,1] or slightly outside it (by 1e-10 or less).
pub const EDGE_CLIP_ERROR_UV_COORD: f64 = 2.25 * DBL_EPSILON;

// faceClipErrorUVCoord is the maximum angle between a returned vertex
// and the nearest point on the exact edge AB expressed as the maximum error
// in an individual u- or v-coordinate. In other words, for each
// returned vertex there is a point on the exact edge AB whose u- and
// v-coordinates differ from the vertex by at most this amount.
pub const FACE_CLIP_ERROR_UV_COORD: f64 = 9.0 * (1.0 / std::f64::consts::SQRT_2) * DBL_EPSILON;

// CellRelation describes the possible relationships between a target cell
// and the cells of the ShapeIndex. If the target is an index cell or is
// contained by an index cell, it is Indexed. If the target is subdivided
// into one or more index cells, it is Subdivided. Otherwise it is Disjoint.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum CellRelation {
    Indexed,
    Subdivided,
    Disjoint,
}

// cellPadding defines the total error when clipping an edge which comes
// from two sources:
// (1) Clipping the original spherical edge to a cube face (the face edge).
//     The maximum error in this step is faceClipErrorUVCoord.
// (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
//     The maximum error in this step is edgeClipErrorUVCoord.
// Finally, since we encounter the same errors when clipping query edges, we
// double the total error so that we only need to pad edges during indexing
// and not at query time.
pub const CELL_PADDING: f64 = 2.0 * (FACE_CLIP_ERROR_UV_COORD + EDGE_CLIP_ERROR_UV_COORD);

// cellSizeToLongEdgeRatio defines the cell size relative to the length of an
// edge at which it is first considered to be long. Long edges do not
// contribute toward the decision to subdivide a cell further. For example,
// a value of 2.0 means that the cell must be at least twice the size of the
// edge in order for that edge to be counted. There are two reasons for not
// counting long edges: (1) such edges typically need to be propagated to
// several children, which increases time and memory costs without much benefit,
// and (2) in pathological cases, many long edges close together could force
// subdivision to continue all the way to the leaf cell level.
pub const CELL_SIZE_TO_LONG_EDGE_RATIO: f64 = 1.0;

// clippedShape represents the part of a shape that intersects a Cell.
// It consists of the set of edge IDs that intersect that cell and a boolean
// indicating whether the center of the cell is inside the shape (for shapes
// that have an interior).
//
// Note that the edges themselves are not clipped; we always use the original
// edges for intersection tests so that the results will be the same as the
// original shape.
#[derive(Debug, Clone)]
pub struct ClippedShape {
    // shapeID is the index of the shape this clipped shape is a part of.
    pub(crate) shape_id: i32,
    // containsCenter indicates if the center of the CellID this shape has been
    // clipped to falls inside this shape. This is false for shapes that do not
    // have an interior.
    pub(crate) contains_center: bool,

    // edges is the ordered set of ShapeIndex original edge IDs. Edges
    // are stored in increasing order of edge ID.
    pub(crate) edges: Vec<i32>,
}

impl ClippedShape {
    // new returns a new clipped shape for the given shapeID and number of expected edges.
    pub fn new(id: i32, num_edges: usize) -> ClippedShape {
        ClippedShape {
            shape_id: id,
            contains_center: false,
            edges: Vec::with_capacity(num_edges),
        }
    }

    // num_edges returns the number of edges that intersect the CellID of the Cell this was clipped to.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    // contains_edge reports if this clipped shape contains the given edge ID.
    pub fn contains_edge(&self, id: i32) -> bool {
        // Linear search is fast because the number of edges per shape is typically
        // very small (less than 10).
        self.edges.iter().any(|&e| e == id)
    }
}

// ShapeIndexCell stores the index contents for a particular CellID.
#[derive(Debug, Clone)]
pub struct ShapeIndexCell {
    pub(crate) shapes: Vec<ClippedShape>,
}

impl ShapeIndexCell {
    // new creates a new cell that is sized to hold the given number of shapes.
    pub fn new(num_shapes: usize) -> ShapeIndexCell {
        ShapeIndexCell {
            shapes: Vec::with_capacity(num_shapes),
        }
    }

    // num_edges reports the total number of edges in all clipped shapes in this cell.
    pub fn num_edges(&self) -> usize {
        self.shapes.iter().map(|cs| cs.num_edges()).sum()
    }

    // add adds the given clipped shape to this index cell.
    pub fn add(&mut self, c: ClippedShape) {
        // C++ uses a set, so it's ordered and unique. We don't currently catch
        // the case when a duplicate value is added.
        self.shapes.push(c);
    }

    // find_by_shape_id returns the clipped shape that contains the given shapeID,
    // or None if none of the clipped shapes contain it.
    pub fn find_by_shape_id(&self, shape_id: i32) -> Option<&ClippedShape> {
        // Linear search is fine because the number of shapes per cell is typically
        // very small (most often 1), and is large only for pathological inputs
        // (e.g. very deeply nested loops).
        self.shapes
            .iter()
            .find(|clipped| clipped.shape_id == shape_id)
    }
}

// faceEdge and ClippedEdge store temporary edge data while the index is being
// updated.
//
// While it would be possible to combine all the edge information into one
// structure, there are two good reasons for separating it:
//
//  - Memory usage. Separating the two means that we only need to
//    store one copy of the per-face data no matter how many times an edge is
//    subdivided, and it also lets us delay computing bounding boxes until
//    they are needed for processing each face (when the dataset spans
//    multiple faces).
//
//  - Performance. UpdateEdges is significantly faster on large polygons when
//    the data is separated, because it often only needs to access the data in
//    ClippedEdge and this data is cached more successfully.

// faceEdge represents an edge that has been projected onto a given face,
#[derive(Debug, Clone)]
pub struct FaceEdge {
    shape_id: i32,              // The ID of shape that this edge belongs to
    edge_id: i32,               // Edge ID within that shape
    max_level: i32,             // Not desirable to subdivide this edge beyond this level
    has_interior: bool,         // Belongs to a shape that has a dimension of 2
    a: crate::r2::point::Point, // The edge endpoints, clipped to a given face
    b: crate::r2::point::Point,
    edge: Edge, // The original edge.
}

// ClippedEdge represents the portion of that edge that has been clipped to a given Cell.
#[derive(Debug, Clone)]
pub struct ClippedEdge {
    face_edge: FaceEdge, // The original unclipped edge
    bound: R2Rect,       // Bounding box for the clipped portion
}

// ShapeIndexIteratorPos defines the set of possible iterator starting positions. By
// default iterators are unpositioned, since this avoids an extra seek in this
// situation where one of the seek methods (such as Locate) is immediately called.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeIndexIteratorPos {
    // Begin specifies the iterator should be positioned at the beginning of the index.
    Begin,
    // End specifies the iterator should be positioned at the end of the index.
    End,
}

// `[ShapeIndexIterator]` is an iterator that provides low-level access to
// the cells of the index. Cells are returned in increasing order of CellID.
//
//	for it := index.Iterator(); !it.Done(); it.Next() {
//	  fmt.Print(it.CellID())
//	}
#[derive(Clone)]
pub struct ShapeIndexIterator<'a> {
    index: &'a ShapeIndex,
    position: usize,
    id: CellID,
    cell: Option<&'a ShapeIndexCell>,
}

impl<'a> Iterator for ShapeIndexIterator<'a> {
    type Item = &'a ShapeIndexCell;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done() {
            return None;
        }

        let cell = self.cell;
        // Advance the iterator for the next call
        self.position += 1;
        self.refresh();

        cell
    }
}

impl<'a> ShapeIndexIterator<'a> {
    // new creates a new iterator for the given index. If a starting
    // position is specified, the iterator is positioned at the given spot.
    pub fn new(index: &'a ShapeIndex, pos: Option<ShapeIndexIteratorPos>) -> Self {
        let mut s = ShapeIndexIterator {
            index,
            position: 0,
            id: CellID::sentinel(),
            cell: None,
        };

        if let Some(pos) = pos {
            match pos {
                ShapeIndexIteratorPos::Begin => s.begin(),
                ShapeIndexIteratorPos::End => s.end(),
            }
        }

        s
    }

    // clone returns a copy of this iterator.
    pub fn clone(&self) -> ShapeIndexIterator<'a> {
        ShapeIndexIterator {
            index: self.index,
            position: self.position,
            id: self.id,
            cell: self.cell,
        }
    }

    // cell_id returns the CellID of the current index cell.
    // If s.Done() is true, a value larger than any valid CellID is returned.
    pub fn cell_id(&self) -> CellID {
        self.id
    }

    // index_cell returns the current index cell.
    pub fn index_cell(&self) -> Option<&'a ShapeIndexCell> {
        self.cell
    }

    // center returns the Point at the center of the current position of the iterator.
    pub fn center(&self) -> Point {
        self.cell_id().center_point()
    }

    // begin positions the iterator at the beginning of the index.
    pub fn begin(&mut self) {
        if !self.index.is_fresh() {
            self.index.maybe_apply_updates();
        }
        self.position = 0;
        self.refresh();
    }

    // next positions the iterator at the next index cell.
    pub fn next(&mut self) {
        self.position += 1;
        self.refresh();
    }

    /// `[prev]` advances the iterator to the previous cell in the index and returns true to
    /// indicate it was not yet at the beginning of the index. If the iterator is at the
    /// first cell the call does nothing and returns false.
    pub fn prev(&mut self) -> bool {
        if self.position <= 0 {
            return false;
        }

        self.position -= 1;
        self.refresh();
        true
    }

    // end positions the iterator at the end of the index.
    pub fn end(&mut self) {
        {
            // Get a lock on the index data to get the length
            let data = self.index.index_data.read().unwrap();
            self.position = data.cells.len();
        }

        // Use refresh to update the iterator state, which handles lifetimes correctly
        self.refresh();
    }

    // done reports if the iterator is positioned at or after the last index cell.
    pub fn done(&self) -> bool {
        self.id == CellID::sentinel()
    }

    // refresh updates the stored internal iterator values.
    fn refresh(&mut self) {
        // Get a copy of the current position
        let pos = self.position;
        let index = self.index;

        // Lock the index data
        let data_guard = index.index_data.read().unwrap();

        if pos < data_guard.cells.len() {
            // Get the cell ID at this position
            self.id = data_guard.cells[pos];

            // For cells, we'll need to get a reference with the right lifetime
            // Clone the cell if it exists to avoid lifetime issues
            if let Some(_cell) = data_guard.cell_map.get(&self.id) {
                // Release the lock before setting the cell reference
                drop(data_guard);

                // Re-acquire the lock and get a properly lifetime'd reference
                self.cell = self.index.get_cell_ref(self.id);
            } else {
                self.cell = None;
            }
        } else {
            self.id = CellID::sentinel();
            self.cell = None;
        }
    }

    // seek positions the iterator at the first cell whose ID >= target, or at the
    // end of the index if no such cell exists.
    pub fn seek(&mut self, target: CellID) {
        // Get a lock on the index data to access cells
        let data = self.index.index_data.read().unwrap();

        // Using binary search to find the position
        self.position = data.cells.binary_search(&target).unwrap_or_else(|pos| pos);

        // Release the lock
        drop(data);

        // Use refresh to update the iterator state, which handles lifetimes correctly
        self.refresh();
    }

    // locate_point positions the iterator at the cell that contains the given Point.
    // If no such cell exists, the iterator position is unspecified, and false is returned.
    // The cell at the matched position is guaranteed to contain all edges that might
    // intersect the line segment between target and the cell's center.
    pub fn locate_point(&mut self, p: Point) -> bool {
        // Let I = cellMap.LowerBound(T), where T is the leaf cell containing
        // point P. Then if T is contained by an index cell, then the
        // containing cell is either I or I'. We test for containment by comparing
        // the ranges of leaf cells spanned by T, I, and I'.
        let target = p.into();
        self.seek(target);
        if !self.done() && self.cell_id().range_min() <= target {
            return true;
        }

        if self.prev() && self.cell_id().range_max() >= target {
            return true;
        }
        false
    }

    // locate_cell_id attempts to position the iterator at the first matching index cell
    // in the index that has some relation to the given CellID. Let T be the target CellID.
    // If T is contained by (or equal to) some index cell I, then the iterator is positioned
    // at I and returns Indexed. Otherwise if T contains one or more (smaller) index cells,
    // then the iterator is positioned at the first such cell I and return Subdivided.
    // Otherwise Disjoint is returned and the iterator position is undefined.
    pub fn locate_cell_id(&mut self, target: CellID) -> CellRelation {
        // Let T be the target, let I = cellMap.LowerBound(T.RangeMin()), and
        // let I' be the predecessor of I. If T contains any index cells, then T
        // contains I. Similarly, if T is contained by an index cell, then the
        // containing cell is either I or I'. We test for containment by comparing
        // the ranges of leaf cells spanned by T, I, and I'.
        self.seek(target.range_min());
        if !self.done() {
            if self.cell_id() >= target && self.cell_id().range_min() <= target {
                return CellRelation::Indexed;
            }
            if self.cell_id() <= target.range_max() {
                return CellRelation::Subdivided;
            }
        }
        if self.prev() && self.cell_id().range_max() >= target {
            return CellRelation::Indexed;
        }
        CellRelation::Disjoint
    }
}

// tracker keeps track of which shapes in a given set contain a particular point
// (the focus). It provides an efficient way to move the focus from one point
// to another and incrementally update the set of shapes which contain it. We use
// this to compute which shapes contain the center of every CellID in the index,
// by advancing the focus from one cell center to the next.
//
// Initially the focus is at the start of the CellID space-filling curve. We then
// visit all the cells that are being added to the ShapeIndex in increasing order
// of CellID. For each cell, we draw two edges: one from the entry vertex to the
// center, and another from the center to the exit vertex (where entry and exit
// refer to the points where the space-filling curve enters and exits the cell).
// By counting edge crossings we can incrementally compute which shapes contain
// the cell center. Note that the same set of shapes will always contain the exit
// point of one cell and the entry point of the next cell in the index, because
// either (a) these two points are actually the same, or (b) the intervening
// cells in CellID order are all empty, and therefore there are no edge crossings
// if we follow this path from one cell to the other.
#[derive(Debug, Clone)]
pub struct Tracker {
    is_active: bool,
    a: Point,
    b: Point,
    next_cell_id: CellID,
    crosser: Option<EdgeCrosser>,
    shape_ids: Vec<i32>,

    // Shape ids saved by save_and_clear_state_before. The state is never saved
    // recursively so we don't need to worry about maintaining a stack.
    saved_ids: Vec<i32>,
}

impl Tracker {
    // new returns a new tracker with the appropriate defaults.
    pub fn new() -> Self {
        // As shapes are added, we compute which ones contain the start of the
        // CellID space-filling curve by drawing an edge from OriginPoint to this
        // point and counting how many shape edges cross this edge.
        let mut t = Tracker {
            is_active: false,
            a: Default::default(),
            b: Self::tracker_origin(),
            next_cell_id: CellID::from_face(0).child_begin_at_level(crate::cellid::MAX_LEVEL),
            crosser: None,
            shape_ids: Vec::new(),
            saved_ids: Vec::new(),
        };

        // CellID curve start
        t.draw_to(face_uv_to_xyz(0, -1.0, -1.0).normalize().into());

        t
    }

    // tracker_origin returns the initial focus point when the tracker is created
    // (corresponding to the start of the CellID space-filling curve).
    fn tracker_origin() -> Point {
        // The start of the S2CellId space-filling curve.
        face_uv_to_xyz(0, -1.0, -1.0).normalize().into()
    }

    // focus returns the current focus point of the tracker.
    pub fn focus(&self) -> Point {
        self.b
    }

    // add_shape adds a shape whose interior should be tracked. contains_focus indicates
    // whether the current focus point is inside the shape. Alternatively, if
    // the focus point is in the process of being moved (via moveTo/drawTo), you
    // can also specify contains_focus at the old focus point and call testEdge
    // for every edge of the shape that might cross the current drawTo line.
    // This updates the state to correspond to the new focus point.
    //
    // This requires shape.has_interior
    pub fn add_shape(&mut self, shape_id: i32, contains_focus: bool) {
        self.is_active = true;
        if contains_focus {
            self.toggle_shape(shape_id);
        }
    }

    // move_to moves the focus of the tracker to the given point. This method should
    // only be used when it is known that there are no edge crossings between the old
    // and new focus locations; otherwise use drawTo.
    pub fn move_to(&mut self, b: Point) {
        self.b = b;
    }

    // draw_to moves the focus of the tracker to the given point. After this method is
    // called, test_edge should be called with all edges that may cross the line
    // segment between the old and new focus locations.
    pub fn draw_to(&mut self, b: Point) {
        self.a = self.b;
        self.b = b;
        // TODO: the edge crosser may need an in-place Init method if this gets expensive
        self.crosser = Some(EdgeCrosser::new(&self.a, &self.b));
    }

    // test_edge checks if the given edge crosses the current edge, and if so, then
    // toggle the state of the given shapeID.
    // This requires shape to have an interior.
    pub fn test_edge(&mut self, shape_id: i32, edge: Edge) {
        if let Some(ref mut crosser) = self.crosser {
            if crosser.edge_or_vertex_crossing(&edge.v0, &edge.v1) {
                self.toggle_shape(shape_id);
            }
        }
    }

    // set_next_cell_id is used to indicate that the last argument to moveTo or drawTo
    // was the entry vertex of the given CellID, i.e. the tracker is positioned at the
    // start of this cell. By using this method together with atCellID, the caller
    // can avoid calling moveTo in cases where the exit vertex of the previous cell
    // is the same as the entry vertex of the current cell.
    pub fn set_next_cell_id(&mut self, next_cell_id: CellID) {
        self.next_cell_id = next_cell_id.range_min();
    }

    // at_cell_id reports if the focus is already at the entry vertex of the given
    // CellID (provided that the caller calls setNextCellID as each cell is processed).
    // TODO: This might not be how rust iterators works
    pub fn at_cell_id(&self, cell_id: CellID) -> bool {
        cell_id.range_min() == self.next_cell_id
    }

    // toggle_shape adds or removes the given shapeID from the set of IDs it is tracking.
    fn toggle_shape(&mut self, shape_id: i32) {
        // Most shapeIDs slices are small, so special case the common steps.

        // If there is nothing here, add it.
        if self.shape_ids.is_empty() {
            self.shape_ids.push(shape_id);
            return;
        }

        // If it's the first element, drop it from the slice.
        if self.shape_ids[0] == shape_id {
            self.shape_ids.remove(0);
            return;
        }

        // Binary search for the insertion position
        match self.shape_ids.binary_search(&shape_id) {
            // Found it, so remove it
            Ok(index) => {
                self.shape_ids.remove(index);
            }
            // Not found, so insert it at the right position
            Err(index) => {
                self.shape_ids.insert(index, shape_id);
            }
        }
    }

    // save_and_clear_state_before makes an internal copy of the state for shape ids below
    // the given limit, and then clear the state for those shapes. This is used during
    // incremental updates to track the state of added and removed shapes separately.
    pub fn save_and_clear_state_before(&mut self, limit_shape_id: i32) {
        let limit = self.lower_bound(limit_shape_id);
        self.saved_ids = self.shape_ids[..limit as usize].to_vec();
        self.shape_ids.drain(..limit as usize);
    }

    // restore_state_before restores the state previously saved by save_and_clear_state_before.
    // This only affects the state for shapeIDs below "limitShapeID".
    pub fn restore_state_before(&mut self, limit_shape_id: i32) {
        let limit = self.lower_bound(limit_shape_id);
        let mut new_ids = self.saved_ids.clone();
        new_ids.extend_from_slice(&self.shape_ids[limit as usize..]);
        self.shape_ids = new_ids;
        self.saved_ids.clear();
    }

    // lower_bound returns the index of the first entry x where x >= shapeID.
    fn lower_bound(&self, shape_id: i32) -> i32 {
        match self.shape_ids.binary_search(&shape_id) {
            Ok(index) => index as i32,
            Err(index) => index as i32,
        }
    }
}

// RemovedShape represents a set of edges from the given shape that is queued for removal.
#[derive(Debug, Clone)]
pub struct RemovedShape {
    shape_id: i32,
    has_interior: bool,
    contains_tracker_origin: bool,
    edges: Vec<Edge>,
}

// Constants for the ShapeIndex status
#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub enum Status {
    Stale,    // There are pending updates.
    Updating, // Updates are currently being applied.
    #[default]
    Fresh, // There are no pending updates.
}

// ShapeIndex indexes a set of Shapes, where a Shape is some collection of edges
// that optionally defines an interior. It can be used to represent a set of
// points, a set of polylines, or a set of polygons. For Shapes that have
// interiors, the index makes it very fast to determine which Shape(s) contain
// a given point or region.
//
// The index can be updated incrementally by adding or removing shapes. It is
// designed to handle up to hundreds of millions of edges. All data structures
// are designed to be small, so the index is compact; generally it is smaller
// than the underlying data being indexed. The index is also fast to construct.
//
// Polygon, Loop, and Polyline implement Shape which allows these objects to
// be indexed easily. You can find useful query methods in CrossingEdgeQuery
// and ClosestEdgeQuery (Not yet implemented in Go).
//
// Example showing how to build an index of Polylines:
//
//	index := NewShapeIndex()
//	for _, polyline := range polylines {
//	    index.Add(polyline);
//	}
//	// Now you can use a CrossingEdgeQuery or ClosestEdgeQuery here.
// #[derive(Debug)]
#[derive(Clone, Default)]
pub struct ShapeIndex {
    // RwLock protected data
    index_data: Arc<RwLock<ShapeIndexData>>,

    // The current status of the index; accessed atomically.
    status: Arc<RwLock<Status>>,

    // The maximum number of edges per cell.
    max_edges_per_cell: i32,
}

impl Debug for ShapeIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{ShapeIndex:{{Status:{:?}, MaxEdgesPerCell: {}}}}}",
            &self.status.read(),
            self.max_edges_per_cell
        )
    }
}

impl ShapeIndex {
    pub(crate) fn num_shapes(&self) -> usize {
        self.index_data.read().unwrap().shapes.len()
    }
}

// ShapeIndexData contains all the data fields for the shape index.
// These fields are protected by a mutex in the ShapeIndex struct.
#[derive(Default, Debug)]
struct ShapeIndexData {
    // shapes is a map of shape ID to shape.
    shapes: HashMap<i32, ShapeType>,

    // nextID tracks the next ID to hand out. IDs are not reused when shapes
    // are removed from the index.
    next_id: i32,

    // cellMap is a map from CellID to the set of clipped shapes that intersect that
    // cell. The cell IDs cover a set of non-overlapping regions on the sphere.
    cell_map: BTreeMap<CellID, ShapeIndexCell>,

    // Track the ordered list of cell IDs.
    cells: Vec<CellID>,

    // pendingAdditionsPos is the index of the first entry that has not been processed
    // via applyUpdatesInternal.
    pending_additions_pos: i32,

    // The set of shapes that have been queued for removal but not processed yet by
    // applyUpdatesInternal.
    pending_removals: Vec<RemovedShape>,
}

impl ShapeIndex {
    // new creates a new ShapeIndex.
    pub fn new() -> Self {
        ShapeIndex {
            index_data: Arc::from(RwLock::new(ShapeIndexData {
                shapes: HashMap::new(),
                cell_map: BTreeMap::new(),
                cells: Vec::new(),
                next_id: 0,
                pending_additions_pos: 0,
                pending_removals: Vec::new(),
            })),
            max_edges_per_cell: 10,
            status: Arc::from(RwLock::new(Fresh)),
        }
    }

    // iterator returns an iterator for this index.
    pub fn iterator(&self) -> ShapeIndexIterator {
        self.maybe_apply_updates();
        ShapeIndexIterator::new(self, Some(ShapeIndexIteratorPos::Begin))
    }

    // begin positions the iterator at the first cell in the index.
    pub fn begin(&self) -> ShapeIndexIterator {
        self.maybe_apply_updates();
        ShapeIndexIterator::new(self, Some(ShapeIndexIteratorPos::Begin))
    }

    // end positions the iterator at the last cell in the index.
    pub fn end(&self) -> ShapeIndexIterator {
        // TODO(roberts): It's possible that updates could happen to the index between
        // the time this is called and the time the iterators position is used and this
        // will be invalid or not the end. For now, things will be undefined if this
        // happens. See about referencing the IsFresh to guard for this in the future.
        self.maybe_apply_updates();
        ShapeIndexIterator::new(self, Some(ShapeIndexIteratorPos::End))
    }

    // region returns a new ShapeIndexRegion for this ShapeIndex.
    pub fn region(&self) -> ShapeIndexRegion {
        ShapeIndexRegion {
            index: self,
            contains_query: ContainsPointQuery::new(self, VertexModel::SemiOpen),
            iter: self.iterator(),
        }
    }

    // len reports the number of Shapes in this index.
    pub fn len(&self) -> usize {
        let data = self.index_data.read().unwrap();
        data.shapes.len()
    }

    // is_empty reports whether this index is empty (has no shapes).
    pub fn is_empty(&self) -> bool {
        let data = self.index_data.read().unwrap();
        data.shapes.is_empty()
    }

    // reset resets the index to its original state.
    pub fn reset(&mut self) {
        let mut data = self.index_data.write().unwrap();
        data.shapes.clear();
        data.next_id = 0;
        data.cell_map.clear();
        data.cells.clear();

        // Reset the status to Fresh
        *self.status.write().unwrap() = Status::Fresh;
    }

    // num_edges returns the number of edges in this index.
    pub fn num_edges(&self) -> usize {
        let data = self.index_data.read().unwrap();
        data.shapes
            .values()
            .map(|shape| shape.num_edges() as usize)
            .sum()
    }

    // num_edges_up_to returns the number of edges in the given index, up to the given
    // limit. If the limit is encountered, the current running total is returned,
    // which may be more than the limit.
    pub fn num_edges_up_to(&self, limit: usize) -> usize {
        let data = self.index_data.read().unwrap();
        let mut num_edges = 0;
        // We choose to iterate over the shapes in order to match the counting
        // up behavior in C++ and for test compatibility instead of using a
        // more idiomatic range over the shape map.
        for i in 0..=data.next_id {
            if let Some(s) = self.shape(i) {
                num_edges += s.num_edges();
                if num_edges >= limit as i64 {
                    break;
                }
            }
        }

        num_edges as usize
    }

    // shape returns the shape with the given ID, or None if the shape has been removed from the index.
    pub fn shape(&self, id: i32) -> Option<ShapeType> {
        let data = self.index_data.read().unwrap();
        // We need to clone the shape to have ownership and avoid lifetime issues
        if let Some(s) = data.shapes.get(&id) {
            // Clone the box to get a new one with the same shape
            // let shape_of_inner = s.clone();
            Some(s.clone())
        } else {
            None
        }
    }

    // id_for_shape returns the id of the given shape in this index, or -1 if it is
    // not in the index.
    //
    // TODO(roberts): Need to figure out an appropriate way to expose this on a Shape.
    // C++ allows a given S2 type (Loop, Polygon, etc) to be part of multiple indexes.
    // By having each type extend S2Shape which has an id element, they all inherit their
    // own id field rather than having to track it themselves.
    pub fn id_for_shape(&self, shape: &ShapeType) -> i32 {
        let data = self.index_data.read().unwrap();
        for (k, v) in &data.shapes {
            if std::ptr::eq(v, shape) {
                return *k;
            }
        }
        -1
    }

    // add adds the given shape to the index and returns the assigned ID.
    pub fn add(&mut self, shape: &ShapeType) -> i32 {
        let mut data = self.index_data.write().unwrap();
        let id = data.next_id;
        data.shapes.insert(id, shape.clone());
        data.next_id += 1;

        // Mark the index as stale
        *self.status.write().unwrap() = Status::Stale;
        id
    }

    // remove removes the given shape from the index.
    pub fn remove(&mut self, shape: &ShapeType) {
        // The index updates itself lazily because it is much more efficient to
        // process additions and removals in batches.
        let id = self.id_for_shape(shape);

        let mut data = self.index_data.write().unwrap();

        // If the shape wasn't found, it's already been removed or was not in the index.
        if !data.shapes.contains_key(&id) {
            return;
        }

        // Remove the shape from the shapes map.
        data.shapes.remove(&id);

        // We are removing a shape that has not yet been added to the index,
        // so there is nothing else to do.
        if id >= data.pending_additions_pos {
            return;
        }

        let num_edges = shape.num_edges();
        let removed = RemovedShape {
            shape_id: id,
            has_interior: shape.dimension() == 2,
            contains_tracker_origin: shape.reference_point().contained,
            edges: (0..num_edges).map(|e| shape.edge(e)).collect(),
        };

        data.pending_removals.push(removed);

        // Mark the index as stale
        *self.status.write().unwrap() = Status::Stale;
    }

    // build triggers the update of the index. Calls to Add and Release are normally
    // queued and processed on the first subsequent query. This has many advantages,
    // the most important of which is that sometimes there *is* no subsequent
    // query, which lets us avoid building the index completely.
    //
    // This method forces any pending updates to be applied immediately.
    pub fn build(&self) {
        self.maybe_apply_updates();
    }

    // is_fresh reports if there are no pending updates that need to be applied.
    // This can be useful to avoid building the index unnecessarily, or for
    // choosing between two different algorithms depending on whether the index
    // is available.
    //
    // The returned index status may be slightly out of date if the index was
    // built in a different thread. This is fine for the intended use (as an
    // efficiency hint), but it should not be used by internal methods.
    pub fn is_fresh(&self) -> bool {
        let status_bind = self.status.read().unwrap();
        *status_bind == Fresh
    }

    // is_first_update reports if this is the first update to the index.
    fn is_first_update(&self) -> bool {
        // Note that it is not sufficient to check whether cellMap is empty, since
        // entries are added to it during the update process.
        let data = self.index_data.read().unwrap();
        data.pending_additions_pos == 0
    }

    // is_shape_being_removed reports if the shape with the given ID is currently slated for removal.
    fn is_shape_being_removed(&self, shape_id: i32) -> bool {
        // All shape ids being removed fall below the index position of shapes being added.
        let data = self.index_data.read().unwrap();
        shape_id < data.pending_additions_pos
    }

    // maybe_apply_updates checks if the index pieces have changed, and if so, applies pending updates.
    fn maybe_apply_updates(&self) {
        // TODO(roberts): To avoid acquiring and releasing the mutex on every
        // query, we should use atomic operations when testing whether the status
        // is fresh and when updating the status to be fresh. This guarantees
        // that any thread that sees a status of fresh will also see the
        // corresponding index updates.
        let status = *self.status.read().unwrap();
        if status != Fresh {
            self.apply_updates_internal();
            // Update the status to Fresh after applying updates
            *self.status.write().unwrap() = Fresh;
        }
    }

    // apply_updates_internal does the actual work of updating the index by applying all
    // pending additions and removals. It does *not* update the indexes status.
    fn apply_updates_internal(&self) {
        // TODO(roberts): Building the index can use up to 20x as much memory per
        // edge as the final index memory size. If this causes issues, add in
        // batched updating to limit the amount of items per batch to a
        // configurable memory footprint overhead.
        let mut t = Tracker::new();

        // allEdges maps a Face to a collection of faceEdges.
        let mut all_edges: Vec<Vec<FaceEdge>> = vec![Vec::new(); 6];

        let mut data = self.index_data.write().unwrap();

        for p in &data.pending_removals {
            self.remove_shape_internal(p, &mut all_edges, &mut t);
        }

        for id in data.pending_additions_pos..data.next_id {
            self.add_shape_internal(id, &mut all_edges, &mut t);
        }

        for face in 0..6 {
            self.update_face_edges(face, &all_edges[face], &mut t);
        }

        data.pending_removals.clear();
        data.pending_additions_pos = data.next_id;
        // It is the caller's responsibility to update the index status.
    }

    // add_shape_internal clips all edges of the given shape to the six cube faces,
    // adds the clipped edges to the set of allEdges, and starts tracking its
    // interior if necessary.
    fn add_shape_internal(&self, shape_id: i32, all_edges: &mut [Vec<FaceEdge>], t: &mut Tracker) {
        let data = self.index_data.read().unwrap();
        let shape = match data.shapes.get(&shape_id) {
            Some(s) => s,
            None => return, // This shape has already been removed.
        };

        let face_edge = FaceEdge {
            shape_id,
            edge_id: 0,
            max_level: 0,
            has_interior: shape.dimension() == 2,
            a: R2Point::new(0.0, 0.0),
            b: R2Point::new(0.0, 0.0),
            edge: Edge::default(),
        };

        if face_edge.has_interior {
            t.add_shape(shape_id, contains_brute_force(shape.clone(), t.focus()));
        }

        let num_edges = shape.num_edges();
        for e in 0..num_edges {
            let edge = shape.edge(e);

            let mut fe = face_edge.clone();
            fe.edge_id = e as i32;
            fe.edge = edge.clone();
            fe.max_level = max_level_for_edge(edge);
            self.add_face_edge(fe, all_edges);
        }
    }

    // add_face_edge adds the given faceEdge into the collection of all edges.
    fn add_face_edge(&self, fe: FaceEdge, all_edges: &mut [Vec<FaceEdge>]) {
        let a_face = face(&fe.edge.v0.0);
        // See if both endpoints are on the same face, and are far enough from
        // the edge of the face that they don't intersect any (padded) adjacent face.
        if a_face == face(&fe.edge.v1.0) {
            let (x, y) = valid_face_xyz_to_uv(a_face, &fe.edge.v0.0);
            let mut fe = fe.clone();
            fe.a = R2Point::new(x, y);

            let (x, y) = valid_face_xyz_to_uv(a_face, &fe.edge.v1.0);
            fe.b = R2Point::new(x, y);

            let max_uv = 1.0 - CELL_PADDING;
            if fe.a.x.abs() <= max_uv
                && fe.a.y.abs() <= max_uv
                && fe.b.x.abs() <= max_uv
                && fe.b.y.abs() <= max_uv
            {
                all_edges[a_face as usize].push(fe);
                return;
            }
        }

        // Otherwise, we simply clip the edge to all six faces.
        for face in 0..6 {
            if let Some((a_clip, b_clip)) = crate::s2::edge_clipping::clip_to_padded_face(
                &fe.edge.v0,
                &fe.edge.v1,
                face,
                CELL_PADDING,
            ) {
                let mut fe = fe.clone();
                fe.a = a_clip;
                fe.b = b_clip;
                all_edges[face as usize].push(fe);
            }
        }
    }

    // update_bound updates the specified endpoint of the given clipped edge and returns the
    // resulting clipped edge.
    fn update_bound(
        &self,
        edge: &ClippedEdge,
        u_end: usize,
        u: f64,
        v_end: usize,
        v: f64,
    ) -> ClippedEdge {
        let mut c = ClippedEdge {
            face_edge: edge.face_edge.clone(),
            bound: R2Rect::empty(),
        };

        if u_end == 0 {
            c.bound = R2Rect::from_points(&[
                R2Point::new(u, if v_end == 0 { v } else { edge.bound.y.lo }),
                R2Point::new(
                    edge.bound.x.hi,
                    if v_end == 0 { edge.bound.y.hi } else { v },
                ),
            ]);
        } else {
            c.bound = R2Rect::from_points(&[
                R2Point::new(
                    edge.bound.x.lo,
                    if v_end == 0 { v } else { edge.bound.y.lo },
                ),
                R2Point::new(u, if v_end == 0 { edge.bound.y.hi } else { v }),
            ]);
        }

        c
    }

    // clip_u_bound clips the given endpoint (lo=0, hi=1) of the u-axis so that
    // it does not extend past the given value of the given edge.
    fn clip_u_bound(&self, edge: &ClippedEdge, u_end: usize, u: f64) -> ClippedEdge {
        // First check whether the edge actually requires any clipping.
        if u_end == 0 {
            if edge.bound.x.lo >= u {
                return edge.clone();
            }
        } else if edge.bound.x.hi <= u {
            return edge.clone();
        }

        // We interpolate the new v-value from the endpoints of the original edge.
        // This has two advantages: (1) we don't need to store the clipped endpoints
        // at all, just their bounding box; and (2) it avoids the accumulation of
        // roundoff errors due to repeated interpolations. The result needs to be
        // clamped to ensure that it is in the appropriate range.
        let e = &edge.face_edge;
        let v = edge
            .bound
            .y
            .clamp_point(crate::s2::edge_clipping::interpolate_f64(
                u, e.a.x, e.b.x, e.a.y, e.b.y,
            ));

        // Determine which endpoint of the v-axis bound to update. If the edge
        // slope is positive we update the same endpoint, otherwise we update the
        // opposite endpoint.
        let v_end = if (u_end == 1) == ((e.a.x > e.b.x) == (e.a.y > e.b.y)) {
            1
        } else {
            0
        };

        self.update_bound(edge, u_end, u, v_end, v)
    }

    // clip_v_bound clips the given endpoint (lo=0, hi=1) of the v-axis so that
    // it does not extend past the given value of the given edge.
    fn clip_v_bound(&self, edge: &ClippedEdge, v_end: usize, v: f64) -> ClippedEdge {
        if v_end == 0 {
            if edge.bound.y.lo >= v {
                return edge.clone();
            }
        } else if edge.bound.y.hi <= v {
            return edge.clone();
        }

        // We interpolate the new u-value from the endpoints of the original edge.
        let e = &edge.face_edge;
        let u = edge
            .bound
            .x
            .clamp_point(interpolate_f64(v, e.a.y, e.b.y, e.a.x, e.b.x));

        // Determine which endpoint of the u-axis bound to update.
        let u_end = if (v_end == 1) == ((e.a.x > e.b.x) == (e.a.y > e.b.y)) {
            1
        } else {
            0
        };

        self.update_bound(edge, u_end, u, v_end, v)
    }

    // clip_v_axis returns the given edge clipped to within the boundaries of the middle
    // interval along the v-axis, and returns the edges for the lower and upper children.
    fn clip_v_axis(
        &self,
        edge: &ClippedEdge,
        middle: r1::interval::Interval,
    ) -> (Option<ClippedEdge>, Option<ClippedEdge>) {
        if edge.bound.y.hi <= middle.lo {
            // Edge is entirely contained in the lower child.
            return (Some(edge.clone()), None);
        } else if edge.bound.y.lo >= middle.hi {
            // Edge is entirely contained in the upper child.
            return (None, Some(edge.clone()));
        }

        // The edge bound spans both children.
        (
            Some(self.clip_v_bound(edge, 1, middle.hi)),
            Some(self.clip_v_bound(edge, 0, middle.lo)),
        )
    }

    // update_face_edges adds or removes the various edges from the index.
    // An edge is added if shapes[id] is not nil, and removed otherwise.
    fn update_face_edges(&self, face: usize, face_edges: &[FaceEdge], t: &mut Tracker) {
        let num_edges = face_edges.len();
        if num_edges == 0 && t.shape_ids.is_empty() {
            return;
        }

        // Create the initial ClippedEdge for each faceEdge. Additional clipped
        // edges are created when edges are split between child cells. We create
        // an array of clipped edges, so that during the recursion we only need
        // to copy pointers in order to propagate an edge to the correct child.
        let mut clipped_edges = Vec::with_capacity(num_edges);
        let mut bound = R2Rect::empty();

        for face_edge in face_edges {
            let clipped = ClippedEdge {
                face_edge: face_edge.clone(),
                bound: R2Rect::from_points(&[face_edge.a, face_edge.b]),
            };

            bound = bound.union(&clipped.bound);
            clipped_edges.push(clipped);
        }

        // Construct the initial face cell containing all the edges, and then update
        // all the edges in the index recursively.
        let face_id = CellID::from_face(face as u64);
        let pcell = crate::s2::padded_cell::PaddedCell::from_cell_id(face_id, CELL_PADDING);

        let disjoint_from_index = self.is_first_update();
        if num_edges > 0 {
            let shrunk_id = self.shrink_to_fit(&pcell, bound);
            if shrunk_id != pcell.id {
                // All the edges are contained by some descendant of the face cell. We
                // can save a lot of work by starting directly with that cell, but if we
                // are in the interior of at least one shape then we need to create
                // index entries for the cells we are skipping over.
                self.skip_cell_range(
                    face_id.range_min(),
                    shrunk_id.range_min(),
                    t,
                    disjoint_from_index,
                );
                let new_pcell =
                    crate::s2::padded_cell::PaddedCell::from_cell_id(shrunk_id, CELL_PADDING);
                self.update_edges(&new_pcell, &clipped_edges, t, disjoint_from_index);
                self.skip_cell_range(
                    shrunk_id.range_max().next(),
                    face_id.range_max().next(),
                    t,
                    disjoint_from_index,
                );
                return;
            }
        }

        // Otherwise (no edges, or no shrinking is possible), subdivide normally.
        self.update_edges(&pcell, &clipped_edges, t, disjoint_from_index);
    }

    // shrink_to_fit shrinks the PaddedCell to fit within the given bounds.
    fn shrink_to_fit(&self, pcell: &crate::s2::padded_cell::PaddedCell, bound: R2Rect) -> CellID {
        let shrunk_id = pcell.shrink_to_fit(&bound);

        if !self.is_first_update() && shrunk_id != pcell.cell_id() {
            // Don't shrink any smaller than the existing index cells, since we need
            // to combine the new edges with those cells.
            let mut iter = self.iterator();
            if iter.locate_cell_id(shrunk_id) == CellRelation::Indexed {
                return iter.cell_id();
            }
        }
        shrunk_id
    }

    // skip_cell_range skips over the cells in the given range, creating index cells if we are
    // currently in the interior of at least one shape.
    fn skip_cell_range(
        &self,
        begin: CellID,
        end: CellID,
        t: &mut Tracker,
        disjoint_from_index: bool,
    ) {
        // If we aren't in the interior of a shape, then skipping over cells is easy.
        if t.shape_ids.is_empty() {
            return;
        }

        // Otherwise generate the list of cell ids that we need to visit, and create
        // an index entry for each one.
        let skipped = crate::s2::cellunion::CellUnion::from_range(begin, end);
        for cell in skipped.0 {
            let empty_edges: Vec<ClippedEdge> = Vec::new();
            let pcell = crate::s2::padded_cell::PaddedCell::from_cell_id(cell, CELL_PADDING);
            self.update_edges(&pcell, &empty_edges, t, disjoint_from_index);
        }
    }

    // update_edges adds or removes the given edges whose bounding boxes intersect a
    // given cell. disjoint_from_index is an optimization hint indicating that cellMap
    // does not contain any entries that overlap the given cell.
    fn update_edges(
        &self,
        pcell: &crate::s2::padded_cell::PaddedCell,
        edges: &[ClippedEdge],
        t: &mut Tracker,
        disjoint_from_index: bool,
    ) {
        // This function is recursive with a maximum recursion depth of 30 (MAX_LEVEL).

        // Incremental updates are handled as follows. All edges being added or
        // removed are combined together in edges, and all shapes with interiors
        // are tracked using tracker. We subdivide recursively as usual until we
        // encounter an existing index cell. At this point we absorb the index
        // cell as follows:
        //
        //   - Edges and shapes that are being removed are deleted from edges and
        //     tracker.
        //   - All remaining edges and shapes from the index cell are added to
        //     edges and tracker.
        //   - Continue subdividing recursively, creating new index cells as needed.
        //   - When the recursion gets back to the cell that was absorbed, we
        //     restore edges and tracker to their previous state.

        let mut index_cell_absorbed = false;
        let mut disjoint_from_index = disjoint_from_index;

        if !disjoint_from_index {
            // There may be existing index cells contained inside pcell. If we
            // encounter such a cell, we need to combine the edges being updated with
            // the existing cell contents by absorbing the cell.
            let mut iter = self.iterator();
            let r = iter.locate_cell_id(pcell.id);
            if r == CellRelation::Disjoint {
                disjoint_from_index = true;
            } else if r == CellRelation::Indexed {
                // Absorb the index cell by transferring its contents to edges and
                // deleting it. We also start tracking the interior of any new shapes.
                self.absorb_index_cell(pcell, &mut iter, edges, t);
                index_cell_absorbed = true;
                disjoint_from_index = true;
            }
        }

        // If there are existing index cells below us, then we need to keep
        // subdividing so that we can merge with those cells. Otherwise,
        // make_index_cell checks if the number of edges is small enough, and creates
        // an index cell if possible (returning true when it does so).
        if !disjoint_from_index || !self.make_index_cell(pcell, edges, t) {
            // Build up vectors of edges to be passed to each child cell.
            // The (i,j) directions are left (i=0), right (i=1), lower (j=0), and upper (j=1).
            let mut child_edges = vec![vec![vec![]; 2]; 2]; // [i][j]

            // Compute the middle of the padded cell, defined as the rectangle in
            // (u,v)-space that belongs to all four (padded) children. By comparing
            // against the four boundaries of middle we can determine which children
            // each edge needs to be propagated to.
            let middle = pcell.middle();

            // Build up a vector of edges to be passed to each child cell.
            // Note that the vast majority of edges are propagated to a single child.
            for edge in edges {
                if edge.bound.x.hi <= middle.x.lo {
                    // Edge is entirely contained in the two left children.
                    let (a, b) = self.clip_v_axis(edge, middle.y);
                    if let Some(edge) = a {
                        child_edges[0][0].push(edge);
                    }
                    if let Some(edge) = b {
                        child_edges[0][1].push(edge);
                    }
                } else if edge.bound.x.lo >= middle.x.hi {
                    // Edge is entirely contained in the two right children.
                    let (a, b) = self.clip_v_axis(edge, middle.y);
                    if let Some(edge) = a {
                        child_edges[1][0].push(edge);
                    }
                    if let Some(edge) = b {
                        child_edges[1][1].push(edge);
                    }
                } else if edge.bound.y.hi <= middle.y.lo {
                    // Edge is entirely contained in the two lower children.
                    if let a = self.clip_u_bound(edge, 1, middle.x.hi) {
                        child_edges[0][0].push(a);
                    }
                    if let b = self.clip_u_bound(edge, 0, middle.x.lo) {
                        child_edges[1][0].push(b);
                    }
                } else if edge.bound.y.lo >= middle.y.hi {
                    // Edge is entirely contained in the two upper children.
                    if let a = self.clip_u_bound(edge, 1, middle.x.hi) {
                        child_edges[0][1].push(a);
                    }
                    if let b = self.clip_u_bound(edge, 0, middle.x.lo) {
                        child_edges[1][1].push(b);
                    }
                } else {
                    // The edge bound spans all four children. The edge
                    // itself intersects either three or four padded children.
                    if let left = self.clip_u_bound(edge, 1, middle.x.hi) {
                        let (a, b) = self.clip_v_axis(&left, middle.y);
                        if let Some(edge) = a {
                            child_edges[0][0].push(edge);
                        }
                        if let Some(edge) = b {
                            child_edges[0][1].push(edge);
                        }
                    }

                    if let right = self.clip_u_bound(edge, 0, middle.x.lo) {
                        let (a, b) = self.clip_v_axis(&right, middle.y);
                        if let Some(edge) = a {
                            child_edges[1][0].push(edge);
                        }
                        if let Some(edge) = b {
                            child_edges[1][1].push(edge);
                        }
                    }
                }
            }

            // Now recursively update the edges in each child. We call the children in
            // increasing order of CellID so that when the index is first constructed,
            // all insertions into cellMap are at the end (which is much faster).
            for pos in 0..4 {
                let (i, j) = pcell.child_ij(pos);
                if !child_edges[i as usize][j as usize].is_empty() || !t.shape_ids.is_empty() {
                    let child_pcell =
                        crate::s2::padded_cell::PaddedCell::from_parent_ij(pcell, i, j);
                    self.update_edges(
                        &child_pcell,
                        &child_edges[i as usize][j as usize],
                        t,
                        disjoint_from_index,
                    );
                }
            }
        }

        if index_cell_absorbed {
            // Restore the state for any edges being removed that we are tracking.
            t.restore_state_before(self.index_data.read().unwrap().pending_additions_pos);
        }
    }

    // make_index_cell builds an indexCell from the given padded cell and set of edges and adds
    // it to the index. If the cell or edges are empty, no cell is added.
    // Returns true if an index cell was created.
    fn make_index_cell(
        &self,
        p: &crate::s2::padded_cell::PaddedCell,
        edges: &[ClippedEdge],
        t: &mut Tracker,
    ) -> bool {
        // If the cell is empty, no index cell is needed. (In most cases this
        // situation is detected before we get to this point, but this can happen
        // when all shapes in a cell are removed.)
        if edges.is_empty() && t.shape_ids.is_empty() {
            return true;
        }

        // Count the number of edges that have not reached their maximum level yet.
        // Return false if there are too many such edges.
        let mut count = 0;
        for ce in edges {
            if (p.level() as i32) < ce.face_edge.max_level {
                count += 1;
            }

            if count > self.max_edges_per_cell {
                return false;
            }
        }

        // Shift the InteriorTracker focus point to the center of the current cell.
        if t.is_active && !edges.is_empty() {
            if !t.at_cell_id(p.id) {
                t.move_to(p.entry_vertex());
            }
            t.draw_to(p.center());
            self.test_all_edges(edges, t);
        }

        // Allocate and fill a new index cell. To get the total number of shapes we
        // need to merge the shapes associated with the intersecting edges together
        // with the shapes that happen to contain the cell center.
        let c_shape_ids = &t.shape_ids;
        let num_shapes = self.count_shapes(edges, c_shape_ids);
        let mut cell = ShapeIndexCell::new(num_shapes);

        // To fill the index cell we merge the two sources of shapes: edge shapes
        // (those that have at least one edge that intersects this cell), and
        // containing shapes (those that contain the cell center). We keep track
        // of the index of the next intersecting edge and the next containing shape
        // as we go along. Both sets of shape ids are already sorted.
        let mut e_next = 0;
        let mut c_next_idx = 0;
        for _i in 0..num_shapes {
            // Advance to next value base + i
            let mut e_shape_id = i32::MAX; // Sentinel
            let mut c_shape_id = i32::MAX; // Sentinel

            if e_next < edges.len() {
                e_shape_id = edges[e_next].face_edge.shape_id;
            }
            if c_next_idx < c_shape_ids.len() {
                c_shape_id = c_shape_ids[c_next_idx];
            }

            let e_begin = e_next;
            if c_shape_id < e_shape_id {
                // The entire cell is in the shape interior.
                let mut clipped = ClippedShape::new(c_shape_id, 0);
                clipped.contains_center = true;
                c_next_idx += 1;
                cell.add(clipped);
            } else {
                // Count the number of edges for this shape and allocate space for them.
                while e_next < edges.len() && edges[e_next].face_edge.shape_id == e_shape_id {
                    e_next += 1;
                }
                let mut clipped = ClippedShape::new(e_shape_id, e_next - e_begin);
                for e in e_begin..e_next {
                    clipped.edges.push(edges[e].face_edge.edge_id);
                }
                if c_shape_id == e_shape_id {
                    clipped.contains_center = true;
                    c_next_idx += 1;
                }
                cell.add(clipped);
            }
        }

        // Add this cell to the map.
        let cell_id = p.id;
        let mut data = self.index_data.write().unwrap();
        data.cell_map.insert(cell_id, cell);
        data.cells.push(cell_id);

        // Shift the tracker focus point to the exit vertex of this cell.
        if t.is_active && !edges.is_empty() {
            t.draw_to(p.exit_vertex());
            self.test_all_edges(edges, t);
            t.set_next_cell_id(p.id.next());
        }

        true
    }

    // test_all_edges calls the trackers testEdge on all edges from shapes that have interiors.
    fn test_all_edges(&self, edges: &[ClippedEdge], t: &mut Tracker) {
        for edge in edges {
            if edge.face_edge.has_interior {
                t.test_edge(edge.face_edge.shape_id, edge.face_edge.edge.clone());
            }
        }
    }

    // count_shapes reports the number of distinct shapes that are either associated with the
    // given edges, or that are currently stored in the InteriorTracker.
    fn count_shapes(&self, edges: &[ClippedEdge], shape_ids: &[i32]) -> usize {
        let mut count = 0;
        let mut last_shape_id = -1;

        // Index of the current element in the shape_ids list.
        let mut shape_id_idx = 0;

        for edge in edges {
            if edge.face_edge.shape_id == last_shape_id {
                continue;
            }

            count += 1;
            last_shape_id = edge.face_edge.shape_id;

            // Skip over any containing shapes up to and including this one,
            // updating count as appropriate.
            while shape_id_idx < shape_ids.len() {
                let clipped_next = shape_ids[shape_id_idx];
                if clipped_next > last_shape_id {
                    break;
                }
                if clipped_next < last_shape_id {
                    count += 1;
                }
                shape_id_idx += 1;
            }
        }

        // Count any remaining containing shapes.
        count += shape_ids.len() - shape_id_idx;
        count
    }

    // absorb_index_cell absorbs an index cell by transferring its contents to edges
    // and/or "tracker", and then delete this cell from the index. If edges includes
    // any edges that are being removed, this method also updates their
    // InteriorTracker state to correspond to the exit vertex of this cell.
    fn absorb_index_cell(
        &self,
        p: &crate::s2::padded_cell::PaddedCell,
        iter: &mut ShapeIndexIterator,
        edges: &[ClippedEdge],
        t: &mut Tracker,
    ) {
        // When we absorb a cell, we erase all the edges that are being removed.
        // However when we are finished with this cell, we want to restore the state
        // of those edges (since that is how we find all the index cells that need
        // to be updated). The edges themselves are restored automatically when
        // UpdateEdges returns from its recursive call, but the InteriorTracker
        // state needs to be restored explicitly.
        //
        // Here we first update the InteriorTracker state for removed edges to
        // correspond to the exit vertex of this cell, and then save the
        // InteriorTracker state. This state will be restored by UpdateEdges when
        // it is finished processing the contents of this cell.
        if t.is_active
            && !edges.is_empty()
            && self.is_shape_being_removed(edges[0].face_edge.shape_id)
        {
            // We probably need to update the tracker. ("Probably" because
            // it's possible that all shapes being removed do not have interiors.)
            if !t.at_cell_id(p.id) {
                t.move_to(p.entry_vertex());
            }
            t.draw_to(p.exit_vertex());
            t.set_next_cell_id(p.id.next());
            for edge in edges {
                let fe = &edge.face_edge;
                if !self.is_shape_being_removed(fe.shape_id) {
                    break; // All shapes being removed come first.
                }
                if fe.has_interior {
                    t.test_edge(fe.shape_id, fe.edge.clone());
                }
            }
        }

        // Save the state of the edges being removed, so that it can be restored
        // when we are finished processing this cell and its children. We don't
        // need to save the state of the edges being added because they aren't being
        // removed from "edges" and will therefore be updated normally as we visit
        // this cell and its children.

        // binding reading

        t.save_and_clear_state_before(self.index_data.read().unwrap().pending_additions_pos);

        // Create a faceEdge for each edge in this cell that isn't being removed.
        let mut face_edges = Vec::new();
        let mut tracker_moved = false;

        if let Some(cell) = iter.index_cell() {
            for clipped in &cell.shapes {
                let shape_id = clipped.shape_id;
                let shape = self.shape(shape_id);
                if shape.is_none() {
                    continue; // This shape is being removed.
                }
                let shape = shape.unwrap();

                let num_clipped = clipped.num_edges();

                // If this shape has an interior, start tracking whether we are inside the
                // shape. updateEdges wants to know whether the entry vertex of this
                // cell is inside the shape, but we only know whether the center of the
                // cell is inside the shape, so we need to test all the edges against the
                // line segment from the cell center to the entry vertex.
                let mut edge = FaceEdge {
                    shape_id,
                    edge_id: 0,
                    max_level: 0,
                    has_interior: shape.dimension() == 2,
                    a: R2Point::new(0.0, 0.0),
                    b: R2Point::new(0.0, 0.0),
                    edge: Edge::default(),
                };

                if edge.has_interior {
                    t.add_shape(shape_id, clipped.contains_center);
                    // There might not be any edges in this entire cell (i.e., it might be
                    // in the interior of all shapes), so we delay updating the tracker
                    // until we see the first edge.
                    if !tracker_moved && num_clipped > 0 {
                        t.move_to(p.center());
                        t.draw_to(p.entry_vertex());
                        t.set_next_cell_id(p.id);
                        tracker_moved = true;
                    }
                }

                for i in 0..num_clipped {
                    let edge_id = clipped.edges[i];
                    edge.edge_id = edge_id;
                    edge.edge = shape.edge(edge_id as i64);
                    edge.max_level = max_level_for_edge(edge.edge.clone());
                    if edge.has_interior {
                        t.test_edge(shape_id, edge.edge.clone());
                    }

                    let clip_result = crate::s2::edge_clipping::clip_to_padded_face(
                        &edge.edge.v0,
                        &edge.edge.v1,
                        p.id.face(),
                        CELL_PADDING,
                    );

                    if let Some((a, b)) = clip_result {
                        edge.a = a;
                        edge.b = b;
                        face_edges.push(edge.clone());
                    }
                }
            }
        }

        // TODO: Complete the implementation of this method

        // Delete this cell from the index.
        let mut data = self.index_data.write().unwrap();
        data.cell_map.remove(&p.id);

        // Remove this cell from the cells vector
        if let Some(pos) = data.cells.iter().position(|&id| id == p.id) {
            data.cells.remove(pos);
        }
    }

    // remove_shape_internal does the actual work for removing a given shape from the index.
    fn remove_shape_internal(
        &self,
        removed: &RemovedShape,
        all_edges: &mut [Vec<FaceEdge>],
        t: &mut Tracker,
    ) {
        let face_edge = FaceEdge {
            shape_id: removed.shape_id,
            edge_id: 0,
            max_level: 0,
            has_interior: removed.has_interior,
            a: R2Point::new(0.0, 0.0),
            b: R2Point::new(0.0, 0.0),
            edge: Edge::default(),
        };

        if face_edge.has_interior {
            t.add_shape(removed.shape_id, removed.contains_tracker_origin);
        }

        for (e_id, edge) in removed.edges.iter().enumerate() {
            let mut fe = face_edge.clone();
            fe.edge_id = e_id as i32;
            fe.edge = edge.clone();
            fe.max_level = max_level_for_edge(fe.edge.clone());
            self.add_face_edge(fe, all_edges);
        }
    }

    // Implement other helper methods as needed for the shape index implementation

    // get_cell_ref returns a properly lifetime'd reference to a cell given its ID
    fn get_cell_ref<'a>(&self, id: CellID) -> Option<&'a ShapeIndexCell> {
        let guard = self.index_data.read().unwrap();
        // Clone the cell if it exists
        if let Some(cell) = guard.cell_map.get(&id) {
            // UNSAFE AS FUCK: We're extending the lifetime from guard to 'a
            // this is only valid if ShapeIndexData lives as long as ShapeIndex
            Some(unsafe { std::mem::transmute::<&ShapeIndexCell, &'a ShapeIndexCell>(cell) })
        } else {
            None
        }
    }
}

// max_level_for_edge reports the maximum level for a given edge.
fn max_level_for_edge(edge: Edge) -> i32 {
    // Compute the maximum cell size for which this edge is considered long.
    // The calculation does not need to be perfectly accurate, so we use Norm
    // rather than Angle for speed.
    let cell_size = edge.v0.sub(Point::from(edge.v1.0)).norm() * CELL_SIZE_TO_LONG_EDGE_RATIO;
    // Now return the first level encountered during subdivision where the
    // average cell size is at most cellSize.
    AVG_EDGEMETRIC.min_level(cell_size) as i32
}

// contains_brute_force determines if the shape contains the given point using a brute force approach.
// This is used for point-in-polygon testing when initializing the Tracker.
fn contains_brute_force(shape: ShapeType, _point: Point) -> bool {
    // If the shape doesn't have an interior, it cannot contain any points
    if shape.dimension() != 2 {
        return false;
    }

    todo!("Contains brute force for Point in ShapeType is NOT implemented")
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_debug_shapeindex() {
        format!("{:?}", ShapeIndex::new());
    }
}
