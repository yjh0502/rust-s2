use crate::cellid::{size_ij, st_to_ij};
use crate::r2::point::Point as R2Point;
use crate::r2::rect::Rect;
use crate::s1::interval::Interval;
use crate::s2;
use crate::s2::cellid::CellID;
use crate::s2::point::Point;
use crate::s2::stuv::{siti_to_st, st_to_uv, uv_to_st};

// Constants used by the PaddedCell.
const SWAP_MASK: u8 = 0x01;
const INVERT_MASK: u8 = 0x02;

// Table to help get the i,j coordinates of the child cell at a given position.
// For each child position pos (0..3) and each parent orientation (0..7),
// posToIJ[orientation][pos] gives the IJNormPos of the child at position pos.
// Each IJNormPos is a 2-bit value where the high bit gives the i-position of
// the child (0=left, 1=right) and the low bit gives the j-position
// (0=bottom, 1=top).
static POS_TO_IJ: [[u8; 4]; 8] = [
    // 0: Normal order: (0,0), (0,1), (1,0), (1,1)
    [0, 1, 2, 3],
    // 1: axes swapped: (0,0), (1,0), (0,1), (1,1)
    [0, 2, 1, 3],
    // 2: i-axis inverted: (1,0), (1,1), (0,0), (0,1)
    [2, 3, 0, 1],
    // 3: i-axis inverted, axes swapped: (1,0), (0,0), (1,1), (0,1)
    [2, 0, 3, 1],
    // 4: j-axis inverted: (0,1), (0,0), (1,1), (1,0)
    [1, 0, 3, 2],
    // 5: j-axis inverted, axes swapped: (0,1), (1,1), (0,0), (1,0)
    [1, 3, 0, 2],
    // 6: i,j axes inverted: (1,1), (1,0), (0,1), (0,0)
    [3, 2, 1, 0],
    // 7: i,j axes inverted, axes swapped: (1,1), (0,1), (1,0), (0,0)
    [3, 1, 2, 0],
];

// Table to invert posToIJ. For each parent orientation (0..7) and
// child ij-position (0..3), ijToPos[orientation][ij] gives the
// position of the child cell (0..3).
static IJ_TO_POS: [[u8; 4]; 8] = [
    // 0: Normal order: (0,0), (0,1), (1,0), (1,1)
    [0, 1, 2, 3],
    // 1: axes swapped: (0,0), (1,0), (0,1), (1,1)
    [0, 2, 1, 3],
    // 2: i-axis inverted: (1,0), (1,1), (0,0), (0,1)
    [2, 3, 0, 1],
    // 3: i-axis inverted, axes swapped: (1,0), (0,0), (1,1), (0,1)
    [1, 3, 0, 2],
    // 4: j-axis inverted: (0,1), (0,0), (1,1), (1,0)
    [1, 0, 3, 2],
    // 5: j-axis inverted, axes swapped: (0,1), (1,1), (0,0), (1,0)
    [2, 0, 3, 1],
    // 6: i,j axes inverted: (1,1), (1,0), (0,1), (0,0)
    [3, 2, 1, 0],
    // 7: i,j axes inverted, axes swapped: (1,1), (0,1), (1,0), (0,0)
    [3, 1, 2, 0],
];

// For each child position (0..3), posToOrientation[pos] gives the change
// in orientation (0..7) to be added to the parent orientation.
static POS_TO_ORIENTATION: [u8; 4] = [
    // (0,0): No change (0)
    0, // (0,1): Output i = input j and output j = input (1-i) (4)
    4, // (1,0): Output i = input (1-j) and output j = input i (0)
    0, // (1,1): Output i = input (1-i) and output j = input (1-j) (4)
    4,
];

// Given a position, posToIJ[orientation][pos] gives the (i,j) coordinates of the child at that
// position in the given orientation. Both i,j are either 0 or 1.
fn pos_to_ij(orientation: u8, pos: u8) -> (u8, u8) {
    let ij = POS_TO_IJ[orientation as usize][pos as usize];
    ((ij >> 1) & 1, ij & 1)
}

// Given a child's (i,j) coordinates, ijToPos[orientation][ij] returns the position of the child.
fn ij_to_pos(orientation: u8, i: u8, j: u8) -> u8 {
    let ij = (2 * i) + j;
    IJ_TO_POS[orientation as usize][ij as usize]
}

// The change in the orientation of the Hilbert curve from the parent to the child is given by
// posToOrientation[pos], where pos is the position of the child (0-3).
fn pos_to_orientation(pos: u8) -> u8 {
    POS_TO_ORIENTATION[pos as usize]
}

/// PaddedCell represents a Cell whose (u,v)-range has been expanded on
/// all sides by a given amount of "padding". Unlike Cell, its methods and
/// representation are optimized for clipping edges against Cell boundaries
/// to determine which cells are intersected by a given set of edges.
#[derive(Debug, Clone)]
pub struct PaddedCell {
    pub(crate) id: CellID,
    padding: f64,
    bound: Rect,
    middle: Option<Rect>, // A rect in (u, v)-space that belongs to all four children.
    i_lo: i64,            // Minimum i-coordinate of this cell before padding
    j_lo: i64,            // Minimum j-coordinate of this cell before padding
    orientation: u8,      // Hilbert curve orientation of this cell.
    level: u8,
}

impl PaddedCell {
    /// Creates a new padded cell with the given padding.
    pub fn from_cell_id(id: CellID, padding: f64) -> Self {
        // Fast path for constructing a top-level face (the most common case).
        if id.is_face() {
            let limit = padding + 1.0;
            let bound =
                Rect::from_points(&[R2Point::new(-limit, limit), R2Point::new(-limit, limit)]);
            let middle = Rect::from_intervals(
                Interval::new(-padding, padding).into(),
                Interval::new(-padding, padding).into(),
            );

            // panic!("ID IS: {}", id);

            return PaddedCell {
                id,
                padding,
                bound,
                middle: Some(middle),
                i_lo: i64::default(),
                j_lo: i64::default(),
                orientation: id.face() & 1,
                level: u8::default(),
            };
        }

        let (_face, i, j, orientation) = id.face_ij_orientation();
        let level = id.level() as u8;
        let ij_size = size_ij(level as u64);
        let i_lo = i as i64 & -(ij_size as i64);
        let j_lo = j as i64 & -(ij_size as i64);

        // Get the bound in (u,v) coordinate space
        let uv_bound = Self::ij_level_to_bound_uv(i, j, level as i64).expanded_by_margin(padding);

        // panic!("{}", id);

        PaddedCell {
            id,
            padding,
            bound: uv_bound,
            middle: None,
            i_lo,
            j_lo,
            orientation: orientation as u8,
            level,
        }
    }

    /// Constructs the child of parent with the given (i,j) index.
    /// The four child cells have indices of (0,0), (0,1), (1,0), (1,1), where the i and j
    /// indices correspond to increasing u- and v-values respectively.
    pub fn from_parent_ij(parent: &PaddedCell, i: u8, j: u8) -> Self {
        // Compute the position and orientation of the child incrementally from the
        // orientation of the parent.
        let pos = ij_to_pos(parent.orientation, i, j);
        let children = parent.id.children();

        let mut cell = PaddedCell {
            id: children[pos as usize],
            padding: parent.padding,
            bound: parent.bound.clone(),
            middle: None,
            orientation: parent.orientation ^ pos_to_orientation(pos),
            level: parent.level + 1,
            i_lo: 0, // Will be set below
            j_lo: 0, // Will be set below
        };

        let ij_size = size_ij(parent.level as u64);
        cell.i_lo = parent.i_lo + ((i as i64) * ij_size as i64) as i64;
        cell.j_lo = parent.j_lo + ((j as i64) * ij_size as i64) as i64;

        // For each child, one corner of the bound is taken directly from the parent
        // while the diagonally opposite corner is taken from middle().
        let middle = parent.middle();
        if i == 1 {
            cell.bound.x.lo = middle.x.lo;
        } else {
            cell.bound.x.hi = middle.x.hi;
        }
        if j == 1 {
            cell.bound.y.lo = middle.y.lo;
        } else {
            cell.bound.y.hi = middle.y.hi;
        }

        cell
    }

    /// Converts an (i,j) coordinate and cell level to a CellID.
    fn ij_level_to_bound_uv(i: i64, j: i64, level: i64) -> Rect {
        let ij_size = size_ij(level as u64);
        let i_lo = i & -(ij_size as i64);
        let j_lo = j & -(ij_size as i64);
        let i_hi = i_lo + (ij_size as i64);
        let j_hi = j_lo + (ij_size as i64);

        // Compute the bound in (s,t) space and convert to (u,v) coordinates.
        let u_lo = st_to_uv(siti_to_st(2 * i_lo as u64));
        let u_hi = st_to_uv(siti_to_st(2 * i_hi as u64));
        let v_lo = st_to_uv(siti_to_st(2 * j_lo as u64));
        let v_hi = st_to_uv(siti_to_st(2 * j_hi as u64));

        Rect::from_intervals(
            Interval::new(u_lo, u_hi).into(),
            Interval::new(v_lo, v_hi).into(),
        )
    }

    /// Returns the CellID this padded cell represents.
    pub fn cell_id(&self) -> CellID {
        self.id
    }

    /// Returns the amount of padding on this cell.
    pub fn padding(&self) -> f64 {
        self.padding
    }

    /// Returns the level this cell is at.
    pub fn level(&self) -> u8 {
        self.level
    }

    /// Returns the center of this cell.
    pub fn center(&self) -> Point {
        let ij_size = size_ij(self.level as u64);
        let si = (2 * self.i_lo + ij_size as i64) as u64;
        let ti = (2 * self.j_lo + ij_size as i64) as u64;
        crate::s2::stuv::face_siti_to_xyz(self.id.face(), si.into(), ti as u64).normalize()
    }

    /// Returns the rectangle in the middle of this cell that belongs to
    /// all four of its children in (u,v)-space.
    pub fn middle(&self) -> Rect {
        // We compute this field lazily because it is not needed the majority of the
        // time (i.e., for cells where the recursion terminates).
        if let Some(ref middle) = self.middle {
            return middle.clone();
        }

        let ij_size = size_ij(self.level as u64);
        let u = st_to_uv(siti_to_st((2 * self.i_lo + ij_size as i64) as u64));
        let v = st_to_uv(siti_to_st((2 * self.j_lo + ij_size as i64) as u64));

        let middle = Rect::from_intervals(
            Interval::new(u - self.padding, u + self.padding).into(),
            Interval::new(v - self.padding, v + self.padding).into(),
        );

        // Store the computed middle rect for future use
        let mut cell = self.clone();
        cell.middle = Some(middle.clone());

        middle
    }

    /// Returns the bounds for this cell in (u,v)-space including padding.
    pub fn bound(&self) -> Rect {
        self.bound.clone()
    }

    /// Returns the (i,j) coordinates for the child cell at the given traversal
    /// position. The traversal position corresponds to the order in which child
    /// cells are visited by the Hilbert curve.
    pub fn child_ij(&self, pos: u8) -> (u8, u8) {
        pos_to_ij(self.orientation, pos)
    }

    /// Returns the vertex where the space-filling curve enters this cell.
    pub fn entry_vertex(&self) -> Point {
        // The curve enters at the (0,0) vertex unless the axis directions are
        // reversed, in which case it enters at the (1,1) vertex.
        let mut i = self.i_lo;
        let mut j = self.j_lo;
        if self.orientation & INVERT_MASK != 0 {
            let ij_size = size_ij(self.level as u64);
            i += ij_size as i64;
            j += ij_size as i64;
        }
        // Convert to u64 before multiplication to avoid overflow
        s2::stuv::face_siti_to_xyz(self.id.face(), (i as u64) * 2, (j as u64) * 2).normalize()
    }

    /// Returns the vertex where the space-filling curve exits this cell.
    pub fn exit_vertex(&self) -> Point {
        // The curve exits at the (1,0) vertex unless the axes are swapped or
        // inverted but not both, in which case it exits at the (0,1) vertex.
        let mut i = self.i_lo;
        let mut j = self.j_lo;
        let ij_size = size_ij(self.level as u64);
        if self.orientation == 0 || self.orientation == SWAP_MASK + INVERT_MASK {
            i += ij_size as i64;
        } else {
            j += ij_size as i64;
        }
        // Convert to u64 before multiplication to avoid overflow
        s2::stuv::face_siti_to_xyz(self.id.face(), (i as u64) * 2, (j as u64) * 2).normalize()
    }

    /// Returns the smallest CellID that contains all descendants of this
    /// padded cell whose bounds intersect the given rect. For algorithms that use
    /// recursive subdivision to find the cells that intersect a particular object, this
    /// method can be used to skip all of the initial subdivision steps where only
    /// one child needs to be expanded.
    ///
    /// Note that this method is not the same as returning the smallest cell that contains
    /// the intersection of this cell with rect. Because of the padding, even if one child
    /// completely contains rect it is still possible that a neighboring child may also
    /// intersect the given rect.
    ///
    /// The provided Rect must intersect the bounds of this cell.
    pub fn shrink_to_fit(&self, rect: &Rect) -> CellID {
        // Quick rejection test: if rect contains the center of this cell along
        // either axis, then no further shrinking is possible.
        if self.level == 0 {
            // Fast path (most calls to this function start with a face cell).
            if rect.x.contains(0.0) || rect.y.contains(0.0) {
                return self.id;
            }
        }

        let ij_size = size_ij(self.level as u64);
        if rect.x.contains(st_to_uv(siti_to_st(
            (2 * self.i_lo + ij_size as i64) as u64,
        ))) || rect.y.contains(st_to_uv(siti_to_st(
            (2 * self.j_lo + ij_size as i64) as u64,
        ))) {
            return self.id;
        }

        // Otherwise we expand rect by the given padding on all sides and find
        // the range of coordinates that it spans along the i- and j-axes. We then
        // compute the highest bit position at which the min and max coordinates
        // differ. This corresponds to the first cell level at which at least two
        // children intersect rect.

        // Increase the padding to compensate for the error in uv_to_st.
        // (The constant below is a provable upper bound on the additional error.)
        const DBL_EPSILON: f64 = 2.2204460492503131e-16;
        let padded = rect.expanded_by_margin(self.padding + 1.5 * DBL_EPSILON);

        let mut i_min = self.i_lo;
        let mut j_min = self.j_lo;
        let i_xor;
        let j_xor;

        // Calculate minimum i-coordinate
        let padded_i_min = st_to_ij(uv_to_st(padded.x.lo)) as i64;
        if i_min < padded_i_min {
            i_min = padded_i_min;
        }

        // Calculate i_xor (XOR of min and max i-coordinates)
        let i_max_limit = self.i_lo + (ij_size as i64) - 1;
        let padded_i_max = st_to_ij(uv_to_st(padded.x.hi));
        if i_max_limit <= padded_i_max {
            i_xor = i_min ^ i_max_limit;
        } else {
            i_xor = i_min ^ padded_i_max;
        }

        // Calculate minimum j-coordinate
        let padded_j_min = st_to_ij(uv_to_st(padded.y.lo));
        if j_min < padded_j_min {
            j_min = padded_j_min;
        }

        // Calculate j_xor (XOR of min and max j-coordinates)
        let j_max_limit = self.j_lo + (ij_size as i64) - 1;
        let padded_j_max = st_to_ij(uv_to_st(padded.y.hi));
        if j_max_limit <= padded_j_max {
            j_xor = j_min ^ j_max_limit;
        } else {
            j_xor = j_min ^ padded_j_max;
        }

        // Compute the highest bit position where the two i- or j-endpoints differ,
        // and then choose the cell level that includes both of these endpoints.
        let level_msb = ((i_xor | j_xor) << 1) + 1;
        let level = if level_msb == 0 {
            30 // MaxLevel
        } else {
            30 - level_msb.leading_zeros() as u8
        };

        if level <= self.level {
            return self.id;
        }

        CellID::from_face_ij(self.id.face(), i_min, j_min).parent(level as u64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::EPSILON;
    use crate::r2::rect::Rect as R2Rect;
    use crate::s2::cellid::CellID;
    use crate::s2::random;

    use rand::Rng;
    use std::f64;

    #[test]
    fn test_padded_cell_from_cell_id() {
        // Test creating a padded cell from a face cell
        let cell_id = CellID::from_face(1);
        let padding = 0.1;
        let padded_cell = PaddedCell::from_cell_id(cell_id, padding);

        assert_eq!(padded_cell.cell_id(), cell_id);
        assert_eq!(padded_cell.padding(), padding);
        assert_eq!(padded_cell.level(), 0);
        assert_eq!(padded_cell.orientation, 1); // Face 1 has orientation 1

        // Verify the bound is correct for a face cell
        let expected_limit = padding + 1.0;
        println!("{:?}", padded_cell);
        assert_eq!(padded_cell.bound.x.lo, -expected_limit);
        assert_eq!(padded_cell.bound.x.hi, expected_limit);
        assert_eq!(padded_cell.bound.y.lo, -expected_limit);
        assert_eq!(padded_cell.bound.y.hi, expected_limit);
    }

    #[test]
    fn test_padded_cell_middle() {
        // Test the middle calculation
        let cell_id = CellID::from_face(0);
        let padding = 0.1;
        let padded_cell = PaddedCell::from_cell_id(cell_id, padding);

        let middle = padded_cell.middle();

        assert_eq!(middle.x.lo, -padding);
        assert_eq!(middle.x.hi, padding);
        assert_eq!(middle.y.lo, -padding);
        assert_eq!(middle.y.hi, padding);
    }

    // #[test]
    // fn test_padded_cell_from_parent_ij() {
    //     // Test creating a child cell from a parent
    //     let parent = PaddedCell::from_cell_id(CellID::from_face(0), 0.1);
    //     let child = PaddedCell::from_parent_ij(&parent, 1, 0);
    //
    //     assert_eq!(child.level(), 1);
    //     assert!(child.cell_id().is_child_of(parent.cell_id()));
    //
    //     // Check that the orientation is updated correctly
    //     assert_eq!(
    //         child.orientation,
    //         parent.orientation ^ pos_to_orientation(ij_to_pos(parent.orientation, 1, 0))
    //     );
    // }

    #[test]
    fn test_trivial_padded_cell_methods() {
        // Test the PaddedCell methods that have approximate Cell equivalents.
        let cid = CellID(0x85dac80c08257ff0);
        println!("{:#x}", cid.0);
        let padding = (1e-15_f64).powf(0.9);
        let cell = s2::cell::Cell::from(cid);
        let p_cell = PaddedCell::from_cell_id(cid, padding);

        println!("cell_id mismatch? {} != {}", cell.id, p_cell.id);
        assert_eq!(cell.id, p_cell.id, "cell_id mismatch");
        println!(
            "level mismatch? {} != {}",
            cell.id.level(),
            p_cell.level() as u64
        );
        assert_eq!(cell.id.level(), p_cell.level() as u64, "level mismatch");
        println!("padding mismatch? {} != {}", padding, p_cell.padding());
        assert_eq!(padding, p_cell.padding(), "padding mismatch");

        println!(
            "bound mismatch? {:?} != {:?}",
            p_cell.bound(),
            cell.bound_uv().expanded_by_margin(padding)
        );
        assert_eq!(
            p_cell.bound(),
            cell.bound_uv().expanded_by_margin(padding),
            "bound mismatch"
        );

        let r = R2Rect::from_points(&[cell.id.center_uv()]).expanded_by_margin(padding);
        println!("middle mismatch? {:?} != {:?}", r, p_cell.middle());
        assert_eq!(r, p_cell.middle(), "middle mismatch");

        println!(
            "center mismatch? {:?} != {:?}",
            cell.id.center_point(),
            p_cell.center()
        );
        assert_eq!(cell.id.center_point(), p_cell.center(), "center mismatch");

        if cid.is_leaf() {
            return;
        }
        let children = cell.children().unwrap();
        for pos in 0..4 {
            let (i, j) = p_cell.child_ij(pos);
            println!("{}th child now!", pos);

            let cell_child = &children[pos as usize];
            let p_cell_child = PaddedCell::from_parent_ij(&p_cell, i, j);

            assert_eq!(cell_child.id, p_cell_child.id, "child cell_id mismatch");
            assert_eq!(
                cell_child.id.level(),
                p_cell_child.level() as u64,
                "child level mismatch"
            );
            assert_eq!(padding, p_cell_child.padding(), "child padding mismatch");

            println!("child x lo: {:.64}", p_cell_child.bound().x.lo);
            println!("child x.hi: {:.64}", p_cell_child.bound().x.hi);
            println!("child y lo: {:.64}", p_cell_child.bound().y.lo);
            println!("child y hi: {:.64}", p_cell_child.bound().y.hi);

            println!(
                "cell child bound with margin expansion - x lo: {:.64}",
                cell_child.bound_uv().expanded_by_margin(padding).x.lo
            );
            println!(
                "cell child bound with margin expansion - x.hi: {:.64}",
                cell_child.bound_uv().expanded_by_margin(padding).x.hi
            );
            println!(
                "cell child bound with margin expansion - y lo: {:.64}",
                cell_child.bound_uv().expanded_by_margin(padding).y.lo
            );
            println!(
                "cell child bound with margin expansion - y hi: {:.64}",
                cell_child.bound_uv().expanded_by_margin(padding).y.hi
            );

            assert_f64_eq!(
                p_cell_child.bound().x.lo,
                cell_child.bound_uv().expanded_by_margin(padding).x.lo
            );
            assert_f64_eq!(
                p_cell_child.bound().x.hi,
                cell_child.bound_uv().expanded_by_margin(padding).x.hi
            );
            assert_f64_eq!(
                p_cell_child.bound().y.lo,
                cell_child.bound_uv().expanded_by_margin(padding).y.lo
            );
            assert_f64_eq!(
                p_cell_child.bound().y.hi,
                cell_child.bound_uv().expanded_by_margin(padding).y.hi
            );

            assert_eq!(
                p_cell_child.bound(),
                cell_child.bound_uv().expanded_by_margin(padding),
                "child bound mismatch"
            );

            let r = R2Rect::from_points(&[cell_child.id.center_uv()]).expanded_by_margin(padding);
            // Using approx_equal instead of direct equality for floating point
            assert!(r.approx_eq(&p_cell_child.middle()), "child middle mismatch");

            assert_eq!(
                cell_child.id.center_point(),
                p_cell_child.center(),
                "child center mismatch"
            );
        }
    }

    #[test]
    fn test_padded_cell_methods() {
        // Test the PaddedCell methods that have approximate Cell equivalents.
        let mut rng = random::rng();

        for _ in 0..1000 {
            let cid = random::cellid(&mut rng);
            println!("{:#x}", cid.0);
            let padding = (1e-15_f64).powf(rng.gen_range(0.0..1.0));
            let cell = s2::cell::Cell::from(cid);
            let p_cell = PaddedCell::from_cell_id(cid, padding);

            assert_eq!(cell.id, p_cell.id, "cell_id mismatch");
            assert_eq!(cell.id.level(), p_cell.level() as u64, "level mismatch");
            assert_eq!(padding, p_cell.padding(), "padding mismatch");

            assert_eq!(
                p_cell.bound(),
                cell.bound_uv().expanded_by_margin(padding),
                "bound mismatch"
            );

            let r = R2Rect::from_points(&[cell.id.center_uv()]).expanded_by_margin(padding);
            assert_eq!(r, p_cell.middle(), "middle mismatch");

            assert_eq!(cell.id.center_point(), p_cell.center(), "center mismatch");

            if cid.is_leaf() {
                continue;
            }

            let children = cell.children().unwrap();
            for pos in 0..4 {
                let (i, j) = p_cell.child_ij(pos);

                let cell_child = &children[pos as usize];
                let p_cell_child = PaddedCell::from_parent_ij(&p_cell, i, j);

                assert_eq!(cell_child.id, p_cell_child.id, "child cell_id mismatch");
                assert_eq!(
                    cell_child.id.level(),
                    p_cell_child.level() as u64,
                    "child level mismatch"
                );
                assert_eq!(padding, p_cell_child.padding(), "child padding mismatch");

                assert_eq!(
                    p_cell_child.bound(),
                    cell_child.bound_uv().expanded_by_margin(padding),
                    "child bound mismatch"
                );

                let r =
                    R2Rect::from_points(&[cell_child.id.center_uv()]).expanded_by_margin(padding);
                // Using approx_equal instead of direct equality for floating point
                assert!(r.approx_eq(&p_cell_child.middle()), "child middle mismatch");

                assert_eq!(
                    cell_child.id.center_point(),
                    p_cell_child.center(),
                    "child center mismatch"
                );
            }
        }
    }

    #[test]
    fn test_trivial_exit_vertices() {
        let id = CellID(0x3000000000000000);
        let unpadded = PaddedCell::from_cell_id(id, 0.0);
        println!("{:?}", unpadded);
        let padded = PaddedCell::from_cell_id(id, 0.5);
        println!("{:?}", padded);

        // Check that entry/exit vertices do not depend on padding
        assert_eq!(
            unpadded.entry_vertex(),
            padded.entry_vertex(),
            "entry vertex should not depend on padding"
        );

        assert_eq!(
            unpadded.exit_vertex(),
            padded.exit_vertex(),
            "exit vertex should not depend on padding"
        );

        // Check that the exit vertex of one cell is the same as the entry vertex
        // of the immediately following cell. This also tests wrapping from the
        // end to the start of the CellID curve with high probability.
        assert_eq!(
            unpadded.exit_vertex(),
            PaddedCell::from_cell_id(id.next_wrap(), 0.0).entry_vertex(),
            "exit vertex should match next cell's entry vertex"
        );
    }

    #[test]
    fn test_padded_cell_entry_exit_vertices() {
        let mut rng = random::rng();

        for _ in 0..1000 {
            let id = random::cellid(&mut rng);
            println!("{:#x}", id.0);
            let unpadded = PaddedCell::from_cell_id(id, 0.0);
            let padded = PaddedCell::from_cell_id(id, 0.5);

            // Check that entry/exit vertices do not depend on padding
            assert_eq!(
                unpadded.entry_vertex(),
                padded.entry_vertex(),
                "entry vertex should not depend on padding"
            );

            assert_eq!(
                unpadded.exit_vertex(),
                padded.exit_vertex(),
                "exit vertex should not depend on padding"
            );

            // Check that the exit vertex of one cell is the same as the entry vertex
            // of the immediately following cell. This also tests wrapping from the
            // end to the start of the CellID curve with high probability.
            assert_eq!(
                unpadded.exit_vertex(),
                PaddedCell::from_cell_id(id.next_wrap(), 0.0).entry_vertex(),
                "exit vertex should match next cell's entry vertex"
            );

            // Check that the entry vertex of a cell is the same as the entry vertex
            // of its first child, and similarly for the exit vertex.
            if id.is_leaf() {
                continue;
            }

            assert_eq!(
                unpadded.entry_vertex(),
                PaddedCell::from_cell_id(id.children()[0], 0.0).entry_vertex(),
                "entry vertex should match first child's entry vertex"
            );

            assert_eq!(
                unpadded.exit_vertex(),
                PaddedCell::from_cell_id(id.children()[3], 0.0).exit_vertex(),
                "exit vertex should match last child's exit vertex"
            );
        }
    }

    #[test]
    fn test_padded_cell_shrink_to_fit() {
        let mut rng = random::rng();

        for _ in 0..1000 {
            // Start with the desired result and work backwards
            let result = random::cellid(&mut rng);
            let result_uv = result.bound_uv();
            let size_uv = result_uv.size();

            // Find the biggest rectangle that fits in "result" after padding
            let max_padding = 0.5 * f64::min(size_uv.x, size_uv.y);
            let padding = max_padding * rng.gen_range(0.0..1.0);
            let max_rect = result_uv.expanded_by_margin(-padding);

            // Start with a random subset of the maximum rectangle
            let mut a = R2Point::new(
                rng.gen_range(max_rect.x.lo..max_rect.x.hi),
                rng.gen_range(max_rect.y.lo..max_rect.y.hi),
            );

            let mut b = R2Point::new(
                rng.gen_range(max_rect.x.lo..max_rect.x.hi),
                rng.gen_range(max_rect.y.lo..max_rect.y.hi),
            );

            if !result.is_leaf() {
                // If the result is not a leaf cell, we must ensure that no child of
                // result also satisfies the conditions of ShrinkToFit(). We do this
                // by ensuring that rect intersects at least two children of result
                // (after padding).
                let use_y = random::one_in(&mut rng, 2);
                let center = if use_y {
                    result.center_uv().y
                } else {
                    result.center_uv().x
                };

                // Find the range of coordinates that are shared between child cells
                // along that axis.
                let shared = if use_y {
                    max_rect
                        .y
                        .intersection(&(center - padding..center + padding).into())
                } else {
                    max_rect
                        .x
                        .intersection(&(center - padding..center + padding).into())
                };

                let mid = rng.gen_range(shared.lo..shared.hi);

                if use_y {
                    a.y = rng.gen_range(max_rect.y.lo..mid);
                    b.y = rng.gen_range(mid..max_rect.y.hi);
                } else {
                    a.x = rng.gen_range(max_rect.x.lo..mid);
                    b.x = rng.gen_range(mid..max_rect.x.hi);
                }
            }

            let rect = R2Rect::from_points(&[a, b]);

            // Choose an arbitrary ancestor as the PaddedCell
            let initial_id = result.parent(rng.gen_range(0..=result.level()));
            let p_cell = PaddedCell::from_cell_id(initial_id, padding);

            assert_eq!(
                p_cell.shrink_to_fit(&rect),
                result,
                "shrink_to_fit returned incorrect cell"
            );
        }
    }
}
