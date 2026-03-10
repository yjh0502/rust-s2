use std;
use std::cmp::PartialEq;
use std::f64::consts::{FRAC_PI_2, PI};

use crate::consts::*;
use crate::r1;
use crate::r3::vector::Vector;
use crate::s1::*;
use crate::s2::edgeutil;
use crate::s2::latlng::LatLng;

#[derive(Clone, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rect {
    pub lat: r1::interval::Interval, // THIS IS R1 INTERVAL
    pub lng: Interval,               // THIS IS S1 INTERVAL
}

impl std::fmt::Debug for Rect {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[lo{:?}, hi{:?}]", self.lo(), self.hi())
    }
}

pub const VALID_RECT_LAT_RANGE: r1::interval::Interval = r1::interval::Interval {
    lo: -PI / 2.,
    hi: PI / 2.,
};
pub const VALID_RECT_LNG_RANGE: Interval = interval::FULL;

impl Rect {
    pub fn lat_lo(&self) -> Angle {
        Angle::from(Rad(self.lat.lo))
    }

    pub fn lat_hi(&self) -> Angle {
        Angle::from(Rad(self.lat.hi))
    }

    pub fn lng_lo(&self) -> Angle {
        Angle::from(Rad(self.lng.lo))
    }

    pub fn lng_hi(&self) -> Angle {
        Angle::from(Rad(self.lng.hi))
    }

    pub fn empty() -> Rect {
        Rect {
            lat: r1::interval::EMPTY,
            lng: interval::EMPTY,
        }
    }

    pub fn full() -> Rect {
        Rect {
            lat: VALID_RECT_LAT_RANGE,
            lng: VALID_RECT_LNG_RANGE,
        }
    }

    pub fn from_center_size(center: LatLng, size: LatLng) -> Self {
        let half = LatLng {
            lat: size.lat * 0.5,
            lng: size.lng * 0.5,
        };
        Rect::from(center).expanded(&half)
    }

    pub fn from_degrees(lat_lo: f64, lng_lo: f64, lat_hi: f64, lng_hi: f64) -> Self {
        Rect {
            lat: r1::interval::Interval::new(
                Angle::from(Deg(lat_lo)).rad(),
                Angle::from(Deg(lat_hi)).rad(),
            ),
            lng: Interval::new(
                Angle::from(Deg(lng_lo)).rad(),
                Angle::from(Deg(lng_hi)).rad(),
            ),
        }
    }

    // Construct the minimal bounding rectangle containing the two given
    // normalized points.  This is equivalent to starting with an empty
    // rectangle and calling add_point() twice.
    pub fn from_point_pair(p1: &LatLng, p2: &LatLng) -> Self {
        Self {
            lat: r1::interval::Interval::from_point_pair(p1.lat.rad(), p2.lat.rad()),
            lng: Interval::from_point_pair(p1.lng.rad(), p2.lng.rad()),
        }
    }

    pub fn is_valid(&self) -> bool {
        self.lat.lo.abs() <= PI / 2.
            && self.lat.hi <= PI / 2.
            && self.lng.is_valid()
            && self.lat.is_empty() == self.lng.is_empty()
    }

    pub fn is_empty(&self) -> bool {
        self.lat.is_empty()
    }

    pub fn is_full(&self) -> bool {
        self.lat == VALID_RECT_LAT_RANGE && self.lng.is_full()
    }

    pub fn is_point(&self) -> bool {
        self.lat.lo == self.lat.hi && self.lng.lo == self.lng.hi
    }

    /// Returns true if `lng.lo()` > `lng.hi()`, i.e. the rectangle crosses
    /// the 180 degree longitude line.
    pub fn is_inverted(&self) -> bool {
        self.lng.is_inverted()
    }

    /// Returns the k-th vertex of the rectangle (k = 0,1,2,3) in
    /// counter-clockwise order (lower left, lower right, upper right,
    /// upper left).  For convenience, the argument is reduced modulo 4
    /// to the range [0..3].
    pub fn vertex(&self, i: u8) -> LatLng {
        let (lat, lng) = match i % 4 {
            0 => (self.lat.lo, self.lng.lo),
            1 => (self.lat.lo, self.lng.hi),
            2 => (self.lat.hi, self.lng.hi),
            3 => (self.lat.hi, self.lng.lo),
            _ => unimplemented!(),
        };
        LatLng {
            lat: Rad(lat).into(),
            lng: Rad(lng).into(),
        }
    }

    pub fn lo(&self) -> LatLng {
        self.vertex(0)
    }

    pub fn hi(&self) -> LatLng {
        self.vertex(2)
    }

    pub fn center(&self) -> LatLng {
        LatLng {
            lat: Rad(self.lat.center()).into(),
            lng: Rad(self.lng.center()).into(),
        }
    }

    pub fn size(&self) -> LatLng {
        LatLng {
            lat: Rad(self.lat.len()).into(),
            lng: Rad(self.lng.len()).into(),
        }
    }

    pub fn area(&self) -> f64 {
        if self.is_empty() {
            0.
        } else {
            let cap_diff = (self.lat.hi.sin() - self.lat.lo.sin()).abs();
            self.lng.len() * cap_diff
        }
    }

    pub fn expanded(&self, margin: &LatLng) -> Self {
        let lat = self.lat.expanded(margin.lat.rad());
        let lng = self.lng.expanded(margin.lng.rad());

        if lat.is_empty() || lng.is_empty() {
            Self::empty()
        } else {
            Rect {
                lat: lat.intersection(&VALID_RECT_LAT_RANGE),
                lng,
            }
        }
    }

    pub fn polar_closure(&self) -> Self {
        if self.lat.lo == -PI / 2. || self.lat.hi == PI / 2. {
            Rect {
                lat: self.lat,
                lng: interval::FULL,
            }
        } else {
            self.clone()
        }
    }

    pub fn union(&self, other: &Self) -> Self {
        Rect {
            lat: self.lat.union(&other.lat),
            lng: self.lng.union(&other.lng),
        }
    }

    pub fn intersection(&self, other: &Self) -> Self {
        let lat = self.lat.intersection(&other.lat);
        let lng = self.lng.intersection(&other.lng);

        if lat.is_empty() || lng.is_empty() {
            Self::empty()
        } else {
            Rect { lat, lng }
        }
    }

    pub fn intersects(&self, other: &Rect) -> bool {
        self.lat.intersects(&other.lat) && self.lng.intersects(&other.lng)
    }

    /// Returns true if the boundary of this rectangle intersects the given
    /// geodesic edge (v0, v1).
    pub fn boundary_intersects(&self, v0: &Point, v1: &Point) -> bool {
        if self.is_empty() {
            return false;
        }
        if !self.lng.is_full() {
            if intersects_lng_edge(v0, v1, self.lat, Angle::from(Rad(self.lng.lo))) {
                return true;
            }
            if intersects_lng_edge(v0, v1, self.lat, Angle::from(Rad(self.lng.hi))) {
                return true;
            }
        }
        if self.lat.lo != -FRAC_PI_2
            && intersects_lat_edge(v0, v1, Angle::from(Rad(self.lat.lo)), self.lng)
        {
            return true;
        }
        if self.lat.hi != FRAC_PI_2
            && intersects_lat_edge(v0, v1, Angle::from(Rad(self.lat.hi)), self.lng)
        {
            return true;
        }
        false
    }

    pub fn interior_intersects(&self, other: &Rect) -> bool {
        self.lat.interior_intersects(&other.lat) && self.lng.interior_intersects(&other.lng)
    }

    // extra functions
    pub fn approx_eq(&self, other: &Self) -> bool {
        self.lat.approx_eq(&other.lat) && self.lng.approx_eq(&other.lng)
    }

    pub fn approx_eq_by(&self, other: &Self, max_error: &LatLng) -> bool {
        self.lat.approx_eq_by(&other.lat, max_error.lat.rad())
            && self.lng.approx_eq_by(&other.lng, max_error.lng.rad())
    }

    // distance_to_latlng returns the minimum distance (measured along the surface of the sphere)
    // from a given point to the rectangle (both its boundary and its interior).
    // If self is empty, the result is meaningless.
    // The latlng must be valid.
    pub fn distance_to_latlng(&self, ll: &LatLng) -> Angle {
        if self.lng.contains(ll.lng.rad()) {
            return (ll.lat - self.lat.hi)
                .max(self.lat.lo - ll.lat)
                .max(Rad(0.).into());
        }
        let i = Interval::new(self.lng.hi, self.lng.complement_center());
        let mut rect_lng = self.lng.lo;
        if i.contains(ll.lng.rad()) {
            rect_lng = self.lng.hi;
        }
        let lo = LatLng {
            lat: Rad(self.lat.lo).into(),
            lng: Rad(rect_lng).into(),
        };
        let hi = LatLng {
            lat: Rad(self.lat.hi).into(),
            lng: Rad(rect_lng).into(),
        };
        distance_from_segment(&Point::from(ll), &Point::from(lo), &Point::from(hi))
    }

    /// hausdorff_distance returns the undirected Hausdorff distance
    /// (measured along the surface of the sphere) to the given rect.
    /// The Hausdorff distance between rectangle A and rectangle B is given by
    ///     H(A, B) = max{h(A, B), h(B, A)}.
    pub fn hausdorff_distance(&self, other: &Self) -> Angle {
        self.directed_hausdorff_distance(other)
            .max(other.directed_hausdorff_distance(self))
    }

    /// directed_hausdorff_distance returns the directed Hausdorff distance
    /// (measured along the surface of the sphere) to the given Rect.
    /// The directed Hausdorff distance from rectangle A to rectangle B
    /// is given by
    ///     h(A, B) = max_{p in A} min_{q in B} d(p, q).
    pub fn directed_hausdorff_distance(&self, other: &Self) -> Angle {
        if self.is_empty() {
            return Angle::from(Rad(0.));
        }
        if other.is_empty() {
            return Angle::from(Rad(PI));
        }
        let lng_distance = self.lng.directed_hausdorff_distance(&other.lng);
        return Self::hausdorff_distance_helper(lng_distance, &self.lat, &other.lat);
    }

    /// Return the directed Hausdorff distance from one longitudinal edge
    /// spanning latitude range `a` to the other longitudinal edge
    /// spanning latitude range `b`, with their longitudinal difference
    /// given by `lng_diff`.
    fn hausdorff_distance_helper(
        lng_diff: Angle,
        a: &r1::interval::Interval,
        b: &r1::interval::Interval,
    ) -> Angle {
        // By symmetry, we can assume a's longtitude is 0 and b's
        // longtitude is lng_diff. Call b's two endpoints b_lo and
        // b_hi. Let H be the hemisphere containing a and delimited by
        // the longitude line of b. The Voronoi diagram of b on H has
        // three edges (portions of great circles) all orthogonal to b
        // and meeting at b_lo cross b_hi.
        //
        // E1: (b_lo, b_lo cross b_hi)
        // E2: (b_hi, b_lo cross b_hi)
        // E3: (-b_mid, b_lo cross b_hi), where b_mid is the midpoint of b
        //
        // They subdivide H into three Voronoi regions. Depending on how
        // longitude 0 (which contains edge a) intersects these regions,
        // we distinguish two cases:
        //
        // Case 1: it intersects three regions. This occurs when
        //         lng_diff <= M_PI_2.
        // Case 2: it intersects only two regions. This occurs when
        //         lng_diff > M_PI_2.
        //
        // In the first case, the directed Hausdorff distance to edge b
        // can only be realized by the following points on a:
        //
        // A1: two endpoints of a.
        // A2: intersection of a with the equator, if b also intersects
        //     the equator.
        //
        // In the second case, the directed Hausdorff distance to edge b
        // can only be  realized by the following points on a:
        //
        // B1: two endpoints of a.
        // B2: intersection of a with E3.
        // B3: farthest point from b_lo to the interior of D,
        //     and farthest point from b_hi to the interior of U, if any,
        //     where D (resp. U) is the portion of edge a below (resp. above)
        //     the intersection point from B2.
        let lng_diff_rad = lng_diff.rad();
        assert!(lng_diff_rad >= 0.0 && lng_diff_rad <= PI);
        if lng_diff_rad == 0.0 {
            return Angle::from(Rad(a.directed_hausdorff_distance(b)));
        }

        // Assumed longitude of b.
        let b_lng = lng_diff;

        // Two endpoints of b.
        let b_lo = Point::from(LatLng::new(Rad(b.lo).into(), b_lng));
        let b_hi = Point::from(LatLng::new(Rad(b.hi).into(), b_lng));

        // Cases A1 and B1.
        let zero = Angle::default();
        let a_lo = Point::from(LatLng::new(Rad(a.lo).into(), zero));
        let a_hi = Point::from(LatLng::new(Rad(a.hi).into(), zero));
        let mut max_distance = distance_from_segment(&a_lo, &b_lo, &b_hi)
            .max(distance_from_segment(&a_hi, &b_lo, &b_hi));

        if lng_diff_rad <= FRAC_PI_2 {
            // Case A2.
            if a.contains(0.0) && b.contains(0.0) {
                max_distance = max_distance.max(lng_diff);
            }
            return max_distance;
        }

        // Case B2.
        let p = Self::bisector_intersection(b, b_lng.rad());
        let p_lat = p.latitude().rad();
        if a.contains(p_lat) {
            max_distance = max_distance.max(p.distance(&b_lo));
        }

        // Case B3.
        if p_lat > a.lo {
            let interval = r1::interval::Interval::new(a.lo, p_lat.min(a.hi));
            if let Some(dist) = Self::interior_max_distance(&interval, &b_lo) {
                max_distance = max_distance.max(dist);
            }
        }
        if p_lat < a.hi {
            let interval = r1::interval::Interval::new(p_lat.max(a.lo), a.hi);
            if let Some(dist) = Self::interior_max_distance(&interval, &b_hi) {
                max_distance = max_distance.max(dist);
            }
        }

        max_distance
    }

    /// Return the intersection of longitude 0 with the bisector of an edge
    /// on longitude `lng` and spanning latitude range `lat`.
    fn bisector_intersection(lat: &r1::interval::Interval, lng: f64) -> Point {
        let lng = lng.abs();
        let lat_center = lat.center();

        // A vector orthogonal to the bisector of the given longitudinal edge.
        let ortho_bisector = if lat_center >= 0.0 {
            LatLng::new(
                Angle::from(Rad(lat_center - FRAC_PI_2)),
                Angle::from(Rad(lng)),
            )
        } else {
            LatLng::new(
                Angle::from(Rad(-lat_center - FRAC_PI_2)),
                Angle::from(Rad(lng - PI)),
            )
        };

        // A vector orthogonal to longitude 0.
        let ortho_lng = Point::from_coords(0.0, -1.0, 0.0);
        ortho_lng.cross(&Point::from(ortho_bisector))
    }

    /// Return the max distance from a point b to the segment spanning
    /// latitude range a_lat on longitude 0, if the max occurs in the
    /// interior of a_lat. Otherwise return None.
    fn interior_max_distance(a_lat: &r1::interval::Interval, b: &Point) -> Option<Angle> {
        // Longitude 0 is in the y=0 plane. b.x >= 0 implies that
        // the maximum does not occur in the interior of a_lat.
        if a_lat.is_empty() || b.0.x >= 0.0 {
            return None;
        }

        // Project b to the y=0 plane. The antipodal of the normalized
        // projection is the point at which the maxium distance from b
        // occurs, if it is contained in a_lat.
        let intersection_point = Point::from_coords(-b.0.x, 0.0, -b.0.z);
        let intersection_lat = LatLng::from(intersection_point).lat;
        if a_lat.interior_contains(intersection_lat.rad()) {
            return Some(b.distance(&intersection_point));
        } else {
            return None;
        }
    }
}

impl<'a, 'b> std::ops::Add<&'a LatLng> for &'b Rect {
    type Output = Rect;
    fn add(self, ll: &'a LatLng) -> Self::Output {
        if !ll.is_valid() {
            self.clone()
        } else {
            Rect {
                lat: self.lat + ll.lat.rad(),
                lng: self.lng + ll.lng.rad(),
            }
        }
    }
}

impl PartialEq for Rect {
    fn eq(&self, other: &Self) -> bool {
        self.lat == other.lat && self.lng == other.lng
    }
}

impl From<LatLng> for Rect {
    fn from(ll: LatLng) -> Self {
        Self {
            lat: r1::interval::Interval::from_point(ll.lat.rad()),
            lng: Interval {
                lo: ll.lng.rad(),
                hi: ll.lng.rad(),
            },
        }
    }
}

use crate::s2::cap::Cap;
use crate::s2::cell::Cell;
use crate::s2::edgeutil::distance_from_segment;
use crate::s2::point::Point;
use crate::s2::region::Region;

impl Region for Rect {
    /// cap_bound returns a cap that countains Rect.
    fn cap_bound(&self) -> Cap {
        // We consider two possible bounding caps, one whose axis passes
        // through the center of the lat-long rectangle and one whose axis
        // is the north or south pole.  We return the smaller of the two caps.

        if self.is_empty() {
            return Cap::empty();
        }

        let (pole_z, pole_angle) = if self.lat.hi + self.lat.lo < 0. {
            // South pole axis yields smaller cap.
            (-1., PI / 2. + self.lat.hi)
        } else {
            (1., PI / 2. - self.lat.lo)
        };

        let pole_cap =
            Cap::from_center_angle(&Point::from_coords(0., 0., pole_z), &Rad(pole_angle).into());

        // For bounding rectangles that span 180 degrees or less in longitude, the
        // maximum cap size is achieved at one of the rectangle vertices.  For
        // rectangles that are larger than 180 degrees, we punt and always return a
        // bounding cap centered at one of the two poles.
        if remainder(self.lng.hi - self.lng.lo, 2. * PI) >= 0.
            && self.lng.hi - self.lng.lo < 2. * PI
        {
            let mid_cap = Cap::from(&Point::from(self.center()))
                + Point::from(self.lo())
                + Point::from(self.hi());

            if mid_cap.height() < pole_cap.height() {
                return mid_cap;
            }
        }
        pole_cap
    }

    /// rect_bound returns itself.
    fn rect_bound(&self) -> Rect {
        self.clone()
    }

    /// contains_cell reports whether the given Cell is contained by this Rect.
    fn contains_cell(&self, c: &Cell) -> bool {
        // A latitude-longitude rectangle contains a cell if and only if it contains
        // the cell's bounding rectangle. This test is exact from a mathematical
        // point of view, assuming that the bounds returned by Cell.RectBound()
        // are tight. However, note that there can be a loss of precision when
        // converting between representations -- for example, if an s2.Cell is
        // converted to a polygon, the polygon's bounding rectangle may not contain
        // the cell's bounding rectangle. This has some slightly unexpected side
        // effects; for instance, if one creates an s2.Polygon from an s2.Cell, the
        // polygon will contain the cell, but the polygon's bounding box will not.
        self.contains(&c.rect_bound())
    }

    /// intersects_cell reports whether this rectangle intersects the
    /// given cell. This is an exact test and may be fairly expensive.
    fn intersects_cell(&self, cell: &Cell) -> bool {
        // First we eliminate the cases where one region completely
        // contains the other. Once these are disposed of, then the
        // regions will intersect if and only if their boundaries
        // intersect.
        if self.is_empty() {
            return false;
        }

        if self.contains_point(&Point(cell.id.raw_point())) {
            return true;
        }

        if cell.contains_point(&Point::from(self.center())) {
            return true;
        }

        // Quick rejection test (not required for correctness).
        if !self.intersects(&cell.rect_bound()) {
            return false;
        }

        // Precompute the cell vertices as points and latitude-longitudes. We also
        // check whether the Cell contains any corner of the rectangle, or
        // vice-versa, since the edge-crossing tests only check the edge interiors.
        let mut vertices = Vec::with_capacity(4);
        let mut latlngs = Vec::with_capacity(4);

        for i in 0..4 {
            vertices.push(cell.vertex(i));
            latlngs.push(LatLng::from(vertices[i]));

            if self.contains_latlng(&latlngs[i]) {
                return true;
            }
            if cell.contains_point(&Point::from(self.vertex(i as u8))) {
                return true;
            }
        }

        // Now check whether the boundaries intersect. Unfortunately, a
        // latitude-longitude rectangle does not have straight edges: two edges
        // are curved, and at least one of them is concave.
        for i in 0..4 {
            let edge_lng = Interval::new(latlngs[i].lng.rad(), latlngs[(i + 1) & 3].lng.rad());
            if !self.lng.intersects(&edge_lng) {
                continue;
            }

            let a = &vertices[i];
            let b = &vertices[(i + 1) & 3];
            if edge_lng.contains(self.lng.lo)
                && intersects_lng_edge(a, b, self.lat, Rad(self.lng.lo).into())
            {
                return true;
            }
            if edge_lng.contains(self.lng.hi)
                && intersects_lng_edge(a, b, self.lat, Rad(self.lng.hi).into())
            {
                return true;
            }
            if intersects_lat_edge(a, b, Rad(self.lat.lo).into(), self.lng) {
                return true;
            }
            if intersects_lat_edge(a, b, Rad(self.lat.hi).into(), self.lng) {
                return true;
            }
        }

        return false;
    }
}

// intersectsLatEdge reports whether the edge AB intersects the given edge of constant
// latitude. Requires the points to have unit length.
fn intersects_lat_edge(a: &Point, b: &Point, lat: Angle, lng: Interval) -> bool {
    // Unfortunately, lines of constant latitude are curves on
    // the sphere. They can intersect a straight edge in 0, 1, or 2 points.

    // First, compute the normal to the plane AB that points vaguely north.
    let mut z = a.cross(b).normalize();
    if z.0.z < 0. {
        z = Point(z.0 * -1.)
    }

    // Extend this to an orthonormal frame (x,y,z) where x is the direction
    // where the great circle through AB achieves its maximium latitude.
    let y = z.cross(&Point::from_coords(0., 0., 1.)).normalize();
    let x = y.cross(&z).normalize();

    // Compute the angle "theta" from the x-axis (in the x-y plane defined
    // above) where the great circle intersects the given line of latitude.
    let sin_lat = lat.rad().sin();
    if sin_lat.abs() >= x.0.z {
        // The great circle does not reach the given latitude.
        return false;
    }

    let cos_theta = sin_lat / x.0.z;
    let sin_theta = (1. - cos_theta * cos_theta).sqrt();
    let theta = sin_theta.atan2(cos_theta);

    // The candidate intersection points are located +/- theta in the x-y
    // plane. For an intersection to be valid, we need to check that the
    // intersection point is contained in the interior of the edge AB and
    // also that it is contained within the given longitude interval "lng".

    // Compute the range of theta values spanned by the edge AB.
    let ab_theta = Interval::from_point_pair(
        a.0.dot(&y.0).atan2(a.0.dot(&x.0)),
        b.0.dot(&y.0).atan2(b.0.dot(&x.0)),
    );

    if ab_theta.contains(theta) {
        // Check if the intersection point is also in the given lng interval.
        let isect = (x * cos_theta) + (y * sin_theta);
        if lng.contains(isect.0.y.atan2(isect.0.x)) {
            return true;
        }
    }

    if ab_theta.contains(-theta) {
        // Check if the other intersection point is also in the given lng interval.
        let isect = (x * cos_theta) - (y * sin_theta);
        if lng.contains(isect.0.y.atan2(isect.0.x)) {
            return true;
        }
    }

    return false;
}

fn intersects_lng_edge(a: &Point, b: &Point, lat: r1::interval::Interval, lng: Angle) -> bool {
    // The nice thing about edges of constant longitude is that
    // they are straight lines on the sphere (geodesics).
    edgeutil::simple_crossing(
        a,
        b,
        &Point::from(LatLng::new(Rad(lat.lo).into(), lng)),
        &Point::from(LatLng::new(Rad(lat.hi).into(), lng)),
    )
}

impl Rect {
    /// contains reports whether this Rect contains the other Rect.
    pub fn contains(&self, other: &Self) -> bool {
        self.lat.contains_interval(&other.lat) && self.lng.contains_interval(&other.lng)
    }

    /// contains_latlng reports whether the given LatLng is within the Rect.
    pub fn contains_latlng(&self, ll: &LatLng) -> bool {
        ll.is_valid() && self.lat.contains(ll.lat.rad()) && self.lng.contains(ll.lng.rad())
    }

    /// contains_point reports whether the given Point is within the Rect.
    pub fn contains_point(&self, p: &Point) -> bool {
        self.contains_latlng(&LatLng::from(p))
    }

    /// Returns true if and only if the interior of this rectangle
    /// contains all points of the given other rectangle (including
    /// its boundary).
    pub fn interior_contains(&self, other: &Self) -> bool {
        self.lat.interior_contains_interval(&other.lat)
            && self.lng.interior_contains_interval(&other.lng)
    }

    /// Returns true if and only if the given latlng is contained in
    /// the interior of the region (i.e. the region excluding its
    /// boundary). The argument must be normalized.
    pub fn interior_contains_latlng(&self, ll: &LatLng) -> bool {
        self.lat.interior_contains(ll.lat.rad()) && self.lng.interior_contains(ll.lng.rad())
    }

    /// Returns true if and only if the given point is contained in
    /// the interior of the region (i.e. the region excluding its
    /// boundary).  The point `p` does not need to be normalized.
    pub fn interior_contains_point(&self, p: &Point) -> bool {
        self.interior_contains_latlng(&LatLng::from(p))
    }

    // Centroid returns the true centroid of the given Rect multiplied
    // by its surface area. The result is not unit length, so you may
    // want to normalize it.  Note that in general the centroid is
    // *not* at the center of the rectangle, and in fact it may not
    // even be contained by the rectangle. (It is the “center of mass”
    // of the rectangle viewed as subset of the unit sphere, i.e. it is
    // the point in space about which this curved shape would rotate.)
    //
    // The reason for multiplying the result by the rectangle area is
    // to make it easier to compute the centroid of more complicated
    // shapes. The centroid of a union of disjoint regions can be
    // computed simply by adding their centroid results.
    pub fn centroid(&self) -> Point {
        // When a sphere is divided into slices of constant thickness
        // by a set of parallel planes, all slices have the same
        // surface area. This implies that the z-component of the
        // centroid is simply the midpoint of the z-interval spanned
        // by the Rect.
        //
        // Similarly, it is easy to see that the (x,y) of the centroid
        // lies in the plane through the midpoint of the rectangle’s
        // longitude interval.  We only need to determine the distance
        // ”d“ of this point from the z-axis.
        //
        // Let’s restrict our attention to a particular z-value. In
        // this z-plane, the Rect is a circular arc. The centroid of
        // this arc lies on a radial line through the midpoint of
        // the arc, and a distance from the z-axis of
        //
        //     r * (sin(alpha) / alpha)
        //
        // where r = sqrt(1-z^2) is the radius of the arc, and “alpha”
        // is half of the arc length (i.e., the arc covers longitudes
        // [-alpha, alpha]).
        //
        // To find the centroid distance from the z-axis for the entire
        // rectangle, we just need to integrate over the z-interval.
        // This gives
        //
        //    d = Integrate[sqrt(1-z^2)*sin(alpha)/alpha, z1..z2] / (z2 - z1)
        //
        // where [z1, z2] is the range of z-values covered by the rectangle.
        // This simplifies to
        //
        //    d = sin(alpha)/(2*alpha*(z2-z1))*(z2*r2 - z1*r1 + theta2 - theta1)
        //
        // where [theta1, theta2] is the latitude interval, z1=sin(theta1),
        // z2=sin(theta2), r1=cos(theta1), and r2=cos(theta2).
        //
        // Finally, we want to return not the centroid itself, but the
        // centroid scaled by the area of the rectangle. The area of the
        // rectangle is
        //
        //    A = 2 * alpha * (z2 - z1)
        //
        // which fortunately appears in the denominator of “d”.

        if self.is_empty() {
            return Point::default();
        }

        let z1 = self.lat.lo.sin();
        let z2 = self.lat.hi.sin();
        let r1 = self.lat.lo.cos();
        let r2 = self.lat.hi.cos();

        let alpha = 0.5 * self.lng.len();
        let r = alpha.sin() * (r2 * z2 - r1 * z1 + self.lat.len());
        let lng = self.lng.center();
        let z = alpha * (z2 + z1) * (z2 - z1); // scaled by the area

        Point(Vector {
            x: r * lng.cos(),
            y: r * lng.sin(),
            z,
        })
    }
}

impl std::fmt::Display for Rect {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}, {}", self.lo(), self.hi())
    }
}

/*
// TODO: Implement the following, which are present in the C++ version:
//   - AddPoint(const S2LatLng&)
//   - ExpandedByDistance(S1Angle)
//   - GetDistance(const S2LatLngRect&)
//   - MayIntersect(const S2Cell&)
//   - Encode, Decode
*/

#[cfg(test)]
#[allow(non_upper_case_globals)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    use crate::cellid::CellID;
    use crate::predicates::sign;
    use crate::r1;
    use crate::r3::vector::Vector;
    use crate::s2::random;
    use rand::Rng;
    use std::ops::Add;

    #[test]
    fn test_rect_angle_accessors() {
        let r = Rect::from_degrees(1.1, 2.2, 3.3, 4.4);
        assert_f64_eq!(r.lat_lo().deg(), 1.1);
        assert_f64_eq!(r.lng_lo().deg(), 2.2);
        assert_f64_eq!(r.lat_hi().deg(), 3.3);
        assert_f64_eq!(r.lng_hi().deg(), 4.4);
    }

    #[test]
    fn test_rect_empty_and_full() {
        let tests = [
            (&Rect::empty(), true, true, false, false),
            (&Rect::full(), true, false, true, false),
        ];
        for &(rect, valid, empty, full, point) in &tests {
            assert_eq!(rect.is_valid(), valid);
            assert_eq!(rect.is_empty(), empty);
            assert_eq!(rect.is_full(), full);
            assert_eq!(rect.is_point(), point);
        }
    }

    #[test]
    fn test_rect_area() {
        let tests = [
            (&Rect::empty(), 0f64),
            (&Rect::full(), 4.0 * PI),
            (
                &Rect {
                    lat: r1::interval::Interval {
                        lo: 0f64,
                        hi: PI / 2.0,
                    },
                    lng: Interval {
                        lo: 0f64,
                        hi: PI / 2.0,
                    },
                },
                PI / 2.0,
            ),
        ];
        for &(rect, want) in &tests {
            assert_eq!(rect.area(), want);
        }
    }

    #[test]
    fn test_rect_eq() {
        let r = Rect::from_degrees(1., 2., 3., 4.);
        assert!(r == Rect::from_degrees(1., 2., 3., 4.));
        assert!(r != Rect::from_degrees(5., 6., 7., 8.));
    }

    #[test]
    fn test_rect_fmt() {
        assert_eq!(
            format!("{:?}", Rect::full()),
            "[lo[-90.0000000, -180.0000000], hi[90.0000000, 180.0000000]]"
        );
    }

    #[test]
    fn test_rect_from_latlng() {
        let ll = LatLng::from_degrees(23.0, 47.0);
        let got = Rect::from(ll.clone());
        assert!(got.is_point());
        assert_eq!(got.center(), ll);
    }

    #[test]
    fn test_rect_from_center_size() {
        let tests = [(
            LatLng::from_degrees(80.0, 170.0),
            LatLng::from_degrees(40.0, 60.0),
            Rect::from_degrees(60.0, 140.0, 90.0, -160.0),
        )];
        for (center, size, want) in &tests {
            assert!(want.approx_eq(&Rect::from_center_size(center.clone(), size.clone())));
        }
    }

    #[test]
    fn test_rect_add_point() {
        let tests = [
            (
                Rect {
                    lat: r1::interval::EMPTY,
                    lng: interval::EMPTY,
                },
                LatLng::from_degrees(0.0, 0.0),
                Rect::from_degrees(0.0, 0.0, 0.0, 0.0),
            ),
            (
                Rect::from_degrees(0.0, 0.0, 0.0, 0.0),
                LatLng {
                    lat: Rad(0.0).into(),
                    lng: Rad(-PI / 2.0).into(),
                },
                Rect::from_degrees(0.0, -90.0, 0.0, 0.0),
            ),
            (
                Rect::from_degrees(0.0, -90.0, 0.0, 0.0),
                LatLng {
                    lat: Rad(PI / 4.0).into(),
                    lng: Rad(-PI).into(),
                },
                Rect::from_degrees(0.0, -180.0, 45.0, 0.0),
            ),
            (
                Rect::from_degrees(0.0, -180.0, 45.0, 0.0),
                LatLng {
                    lat: Rad(PI / 2.0).into(),
                    lng: Rad(0.0).into(),
                },
                Rect::from_degrees(0.0, -180.0, 90.0, 0.0),
            ),
        ];

        for (input, point, want) in &tests {
            let got = input.add(point);
            assert!(want.approx_eq(&got));
        }
    }

    #[test]
    fn test_rect_from_point_pair() {
        let p1 = LatLng::from_degrees(0.1, 0.2);
        let p2 = LatLng::from_degrees(30.1, 45.5);
        let expected = Rect::empty().add(&p1).add(&p2);
        assert_eq!(Rect::from_point_pair(&p1, &p2), expected);
    }

    #[test]
    fn test_rect_vertex() {
        let r1 = &Rect {
            lat: r1::interval::Interval::new(0.0, PI / 2.0),
            lng: Interval::new(-PI, 0.0),
        };
        let tests = [
            (
                r1,
                0,
                &LatLng {
                    lat: Rad(0.0).into(),
                    lng: Rad(PI).into(),
                },
            ),
            (
                r1,
                1,
                &LatLng {
                    lat: Rad(0.0).into(),
                    lng: Rad(0.0).into(),
                },
            ),
            (
                r1,
                2,
                &LatLng {
                    lat: Rad(PI / 2.0).into(),
                    lng: Rad(0.0).into(),
                },
            ),
            (
                r1,
                3,
                &LatLng {
                    lat: Rad(PI / 2.0).into(),
                    lng: Rad(PI).into(),
                },
            ),
        ];

        for &(r, i, want) in &tests {
            assert_eq!(r.vertex(i), *want);
            assert_eq!(r.vertex(i + 4), *want);
        }
    }

    #[test]
    fn test_rect_vertex_ccw_order() {
        for i in 0..4 {
            let lat = PI / 4.0 * (i as f64 - 2.0);
            let lng = PI / 2.0 * (i as f64 - 2.0) + 0.2;
            let r = Rect {
                lat: r1::interval::Interval {
                    lo: lat,
                    hi: lat + PI / 4.0,
                },
                lng: Interval::new(lng % (2.0 * PI), (lng + PI / 2.0) % (2.0 * PI)),
            };

            for k in 0i8..4 {
                assert!(sign(
                    &Point::from(r.vertex(((k - 1i8) & 3i8) as u8)),
                    &Point::from(r.vertex(k as u8)),
                    &Point::from(r.vertex(((k + 1i8) & 3i8) as u8))
                ));
            }
        }
    }

    #[test]
    fn test_rect_is_inverted() {
        assert_eq!(Rect::from_degrees(0., 22., 7., 88.).is_inverted(), false);
        assert_eq!(Rect::from_degrees(0., 88., 7., 22.).is_inverted(), true);
    }

    #[test]
    fn test_contains_latlng() {
        let tests = [
            (
                Rect::from_degrees(0.0, -180.0, 90.0, 0.0),
                LatLng::from_degrees(30.0, -45.0),
                true,
            ),
            (
                Rect::from_degrees(0.0, -180.0, 90.0, 0.0),
                LatLng::from_degrees(30.0, 45.0),
                false,
            ),
            (
                Rect::from_degrees(0.0, -180.0, 90.0, 0.0),
                LatLng::from_degrees(0.0, -180.0),
                true,
            ),
            (
                Rect::from_degrees(0.0, -180.0, 90.0, 0.0),
                LatLng::from_degrees(90.0, 0.0),
                true,
            ),
        ];

        for (input, ll, want) in &tests {
            assert_eq!(input.contains_latlng(ll), *want);
        }
    }

    #[test]
    fn test_interior_contains_latlng() {
        let eq_m180 = LatLng::from_degrees(0., -180.0);
        let north_pole = LatLng::from_degrees(90., 0.);
        let r = Rect::from(eq_m180).add(&north_pole);
        assert!(r.interior_contains_latlng(&LatLng::from_degrees(30., -45.)));
        assert!(!r.interior_contains_latlng(&LatLng::from_degrees(30., 45.)));
        assert!(!r.interior_contains_latlng(&eq_m180));
        assert!(!r.interior_contains_latlng(&north_pole));
    }

    #[test]
    fn test_interior_contains_point() {
        let r = Rect::from_degrees(0.0, -180.0, 90.0, 0.0);
        assert!(r.interior_contains_point(&Point::from_coords(0.5, -0.3, 0.1)));
        assert!(!r.interior_contains_point(&Point::from_coords(0.5, 0.2, 0.1)));
    }

    /// Helper for test_interval_ops().
    fn verify_interval_ops(
        x: &Rect,
        y: &Rect,
        expected_relations: &str,
        expected_union: &Rect,
        expected_intersection: &Rect,
    ) {
        let mut s = String::with_capacity(4);
        s.push(if x.contains(y) { 'T' } else { 'F' });
        s.push(if x.interior_contains(y) { 'T' } else { 'F' });
        s.push(if x.intersects(y) { 'T' } else { 'F' });
        s.push(if x.interior_intersects(y) { 'T' } else { 'F' });
        assert_eq!(s, expected_relations, "x={:?} y={:?}", x, y);
        assert_eq!(x.union(y), *expected_union, "x={:?} y={:?}", x, y);
        assert_eq!(
            x.intersection(y),
            *expected_intersection,
            "x={:?} y={:?}",
            x,
            y
        );
        if y.size() == LatLng::default() {
            assert_eq!(x.clone().add(&y.lo()), *expected_union);
        }
    }

    /// Tests on Rect::contains(&Rect), Rect::interior_contains(&Rect),
    /// Rect::intersects(&Rect), Rect::interior_intersects(&Rect),
    /// Rect::union(&Rect), and Rect::intersection(&Rect).
    #[test]
    fn test_interval_ops() {
        // Rectangle "r1" covers one-quarter of the sphere.
        let r1 = Rect::from_degrees(0., -180., 90., 0.);

        // Test operations where one rectangle consists of a single point.
        let r1_mid = Rect::from_degrees(45., -90., 45., -90.);
        verify_interval_ops(&r1, &r1_mid, "TTTT", &r1, &r1_mid);

        let req_m180 = Rect::from_degrees(0., -180., 0., -180.);
        verify_interval_ops(&r1, &req_m180, "TFTF", &r1, &req_m180);

        let rnorth_pole = Rect::from_degrees(90., 0., 90., 0.);
        verify_interval_ops(&r1, &rnorth_pole, "TFTF", &r1, &rnorth_pole);

        verify_interval_ops(
            &r1,
            &Rect::from_degrees(-10., -1., 1., 20.),
            "FFTT",
            &Rect::from_degrees(-10., 180., 90., 20.),
            &Rect::from_degrees(0., -1., 1., 0.),
        );
        verify_interval_ops(
            &r1,
            &Rect::from_degrees(-10., -1., 0., 20.),
            "FFTF",
            &Rect::from_degrees(-10., 180., 90., 20.),
            &Rect::from_degrees(0., -1., 0., 0.),
        );
        verify_interval_ops(
            &r1,
            &Rect::from_degrees(-10., 0., 1., 20.),
            "FFTF",
            &Rect::from_degrees(-10., 180., 90., 20.),
            &Rect::from_degrees(0., 0., 1., 0.),
        );
        verify_interval_ops(
            &Rect::from_degrees(-15., -160., -15., -150.),
            &Rect::from_degrees(20., 145., 25., 155.),
            "FFFF",
            &Rect::from_degrees(-15., 145., 25., -150.),
            &Rect::empty(),
        );
        verify_interval_ops(
            &Rect::from_degrees(70., -10., 90., -140.),
            &Rect::from_degrees(60., 175., 80., 5.),
            "FFTT",
            &Rect::from_degrees(60., -180., 90., 180.),
            &Rect::from_degrees(70., 175., 80., 5.),
        );

        // Check that the intersection of two rectangles that overlap
        // in latitude but not longitude is valid, and vice versa.
        verify_interval_ops(
            &Rect::from_degrees(12., 30., 60., 60.),
            &Rect::from_degrees(0., 0., 30., 18.),
            "FFFF",
            &Rect::from_degrees(0., 0., 60., 60.),
            &Rect::empty(),
        );
        verify_interval_ops(
            &Rect::from_degrees(0., 0., 18., 42.),
            &Rect::from_degrees(30., 12., 42., 60.),
            "FFFF",
            &Rect::from_degrees(0., 0., 42., 60.),
            &Rect::empty(),
        );
    }

    #[test]
    fn test_boundary_intersects_empty_rectangle() {
        let rect = Rect::empty();
        let (lo, hi) = (&Point::from(rect.lo()), &Point::from(rect.hi()));
        assert_eq!(rect.boundary_intersects(lo, lo), false);
        assert_eq!(rect.boundary_intersects(lo, hi), false);
    }

    #[test]
    fn test_boundary_intersects_full_rectangle() {
        let rect = Rect::full();
        let (lo, hi) = (&Point::from(rect.lo()), &Point::from(rect.hi()));
        assert_eq!(rect.boundary_intersects(lo, lo), false);
        assert_eq!(rect.boundary_intersects(lo, hi), false);
    }

    fn pt(lat: i32, lng: i32) -> Point {
        Point::from(LatLng::new(
            Angle::from(Deg(lat as f64)),
            Angle::from(Deg(lng as f64)),
        ))
    }

    #[test]
    fn test_boundary_intersects_spherical_lune() {
        // This rectangle only has two non-degenerate sides.
        let rect = Rect::from_degrees(-90., 100., 90., 120.);
        assert_eq!(rect.boundary_intersects(&pt(60, 60), &pt(90, 60)), false);
        assert_eq!(rect.boundary_intersects(&pt(-60, 110), &pt(60, 110)), false);
        assert_eq!(rect.boundary_intersects(&pt(60, 95), &pt(60, 110)), true);
        assert_eq!(rect.boundary_intersects(&pt(60, 115), &pt(80, 125)), true);
    }

    #[test]
    fn test_boundary_intersects_north_hemisphere() {
        // This rectangle only has one non-degenerate side.
        let rect = Rect::from_degrees(0., -180., 90., 180.);
        assert_eq!(
            rect.boundary_intersects(&pt(60, -180), &pt(90, -180)),
            false
        );
        assert_eq!(rect.boundary_intersects(&pt(60, -170), &pt(60, 170)), false);
        assert_eq!(
            rect.boundary_intersects(&pt(-10, -180), &pt(10, -180)),
            true
        );
    }

    #[test]
    fn test_boundary_intersects_south_hemisphere() {
        // This rectangle only has one non-degenerate side.
        let rect = Rect::from_degrees(-90., -180., 0., 180.);
        assert_eq!(
            rect.boundary_intersects(&pt(-90, -180), &pt(-60, -180)),
            false
        );
        assert_eq!(
            rect.boundary_intersects(&pt(-60, -170), &pt(-60, 170)),
            false
        );
        assert_eq!(
            rect.boundary_intersects(&pt(-10, -180), &pt(10, -180)),
            true
        );
    }

    #[test]
    fn test_boundary_intersects_rect_crossing_anti_meridian() {
        let rect = Rect::from_degrees(20., 170., 40., -170.);

        // Check that crossings of all four sides are detected.
        assert_eq!(rect.boundary_intersects(&pt(25, 160), &pt(25, 180)), true);
        assert_eq!(rect.boundary_intersects(&pt(25, -160), &pt(25, -180)), true);
        assert_eq!(rect.boundary_intersects(&pt(15, 175), &pt(30, 175)), true);
        assert_eq!(rect.boundary_intersects(&pt(45, 175), &pt(30, 175)), true);

        // Check that the edges on the opposite side of the sphere but
        // at the same latitude do not intersect the rectangle boundary.
        assert_eq!(rect.boundary_intersects(&pt(25, -20), &pt(25, 0)), false);
        assert_eq!(rect.boundary_intersects(&pt(25, 20), &pt(25, 0)), false);
        assert_eq!(rect.boundary_intersects(&pt(15, -5), &pt(30, -5)), false);
        assert_eq!(rect.boundary_intersects(&pt(45, -5), &pt(30, -5)), false);
    }

    #[test]
    fn test_rect_expanded() {
        let tests = [
            (
                Rect::from_degrees(70.0, 150.0, 80.0, 170.0),
                LatLng::from_degrees(20.0, 30.0),
                Rect::from_degrees(50.0, 120.0, 90.0, -160.0),
            ),
            (
                Rect::empty(),
                LatLng::from_degrees(20.0, 30.0),
                Rect::empty(),
            ),
            (
                Rect::full(),
                LatLng::from_degrees(500.0, 500.0),
                Rect::full(),
            ),
            (
                Rect::from_degrees(-90.0, 170.0, 10.0, 20.0),
                LatLng::from_degrees(30.0, 80.0),
                Rect::from_degrees(-90.0, -180.0, 40.0, 180.0),
            ),
            // Negative margins.
            (
                Rect::from_degrees(10.0, -50.0, 60.0, 70.0),
                LatLng::from_degrees(-10.0, -10.0),
                Rect::from_degrees(20.0, -40.0, 50.0, 60.0),
            ),
            (
                Rect::from_degrees(-20.0, -180.0, 20.0, 180.0),
                LatLng::from_degrees(-10.0, -10.0),
                Rect::from_degrees(-10.0, -180.0, 10.0, 180.0),
            ),
            (
                Rect::from_degrees(-20.0, -180.0, 20.0, 180.0),
                LatLng::from_degrees(-30.0, -30.0),
                Rect::empty(),
            ),
            (
                Rect::from_degrees(-90.0, 10.0, 90.0, 11.0),
                LatLng::from_degrees(-10.0, -10.0),
                Rect::empty(),
            ),
            (
                Rect::from_degrees(-90.0, 10.0, 90.0, 100.0),
                LatLng::from_degrees(-10.0, -10.0),
                Rect::from_degrees(-80.0, 20.0, 80.0, 90.0),
            ),
            (
                Rect::empty(),
                LatLng::from_degrees(-50.0, -500.0),
                Rect::empty(),
            ),
            (
                Rect::full(),
                LatLng::from_degrees(-50.0, -50.0),
                Rect::from_degrees(-40.0, -180.0, 40.0, 180.0),
            ),
            // Mixed margins.
            (
                Rect::from_degrees(10.0, -50.0, 60.0, 70.0),
                LatLng::from_degrees(-10.0, 30.0),
                Rect::from_degrees(20.0, -80.0, 50.0, 100.0),
            ),
            (
                Rect::from_degrees(-20.0, -180.0, 20.0, 180.0),
                LatLng::from_degrees(10.0, -500.0),
                Rect::from_degrees(-30.0, -180.0, 30.0, 180.0),
            ),
            (
                Rect::from_degrees(-90.0, -180.0, 80.0, 180.0),
                LatLng::from_degrees(-30.0, 500.0),
                Rect::from_degrees(-60.0, -180.0, 50.0, 180.0),
            ),
            (
                Rect::from_degrees(-80.0, -100.0, 80.0, 150.0),
                LatLng::from_degrees(30.0, -50.0),
                Rect::from_degrees(-90.0, -50.0, 90.0, 100.0),
            ),
            (
                Rect::from_degrees(0.0, -180.0, 50.0, 180.0),
                LatLng::from_degrees(-30.0, 500.0),
                Rect::empty(),
            ),
            (
                Rect::from_degrees(-80.0, 10.0, 70.0, 20.0),
                LatLng::from_degrees(30.0, -200.0),
                Rect::empty(),
            ),
            (
                Rect::empty(),
                LatLng::from_degrees(100.0, -100.0),
                Rect::empty(),
            ),
            (
                Rect::full(),
                LatLng::from_degrees(100.0, -100.0),
                Rect::full(),
            ),
        ];

        for (input, margin, want) in &tests {
            let got = input.expanded(margin);
            assert!(want.approx_eq(&got));
        }
    }

    #[test]
    fn test_rect_polar_closure() {
        let tests = [
            (
                Rect::from_degrees(-89.0, 0.0, 89.0, 1.0),
                Rect::from_degrees(-89.0, 0.0, 89.0, 1.0),
            ),
            (
                Rect::from_degrees(-90.0, -30.0, -45.0, 100.0),
                Rect::from_degrees(-90.0, -180.0, -45.0, 180.0),
            ),
            (
                Rect::from_degrees(89.0, 145.0, 90.0, 146.0),
                Rect::from_degrees(89.0, -180.0, 90.0, 180.0),
            ),
            (
                Rect::from_degrees(-90.0, -145.0, 90.0, -144.0),
                Rect::full(),
            ),
        ];

        for (r, want) in &tests {
            let got = r.polar_closure();
            assert!(want.approx_eq(&got));
        }
    }

    #[test]
    fn test_rect_cap_bound() {
        let tests = [
            (
                // Bounding cap at center is smaller.
                Rect::from_degrees(-45.0, -45.0, 45.0, 45.0),
                Cap::from_center_height(
                    &Point(Vector {
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    }),
                    0.5,
                ),
            ),
            (
                // Bounding cap at north pole is smaller.
                Rect::from_degrees(88.0, -80.0, 89.0, 80.0),
                Cap::from_center_angle(
                    &Point(Vector {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    }),
                    &Angle::from(Deg(2.0)),
                ),
            ),
            (
                // Longitude span > 180 degrees.
                Rect::from_degrees(-30.0, -150.0, -10.0, 50.0),
                Cap::from_center_angle(
                    &Point(Vector {
                        x: 0.0,
                        y: 0.0,
                        z: -1.0,
                    }),
                    &Angle::from(Deg(80.)),
                ),
            ),
        ];

        for (r, want) in &tests {
            let got = r.cap_bound();
            assert!(want.approx_eq(&got));
        }
    }

    #[test]
    fn test_rect_interval_ops() {
        // Rectangle that covers one-quarter of the sphere.
        let rect = Rect::from_degrees(0., -180., 90., 0.);

        // Test operations where one rectangle consists of a single point.
        let rect_mid = Rect::from_degrees(45., -90., 45., -90.);
        let rect180 = Rect::from_degrees(0., -180., 0., -180.);
        let north_pole = Rect::from_degrees(90., 0., 90., 0.);

        struct Test<'a> {
            rect: &'a Rect,
            other: &'a Rect,
            contains: bool,
            intersects: bool,
            union: &'a Rect,
            intersection: &'a Rect,
        }
        let tests: [Test; 10] = [
            Test {
                rect: &rect,
                other: &rect_mid,
                contains: true,
                intersects: true,
                union: &rect,
                intersection: &rect_mid,
            },
            Test {
                rect: &rect,
                other: &rect180,
                contains: true,
                intersects: true,
                union: &rect,
                intersection: &rect180,
            },
            Test {
                rect: &rect,
                other: &north_pole,
                contains: true,
                intersects: true,
                union: &rect,
                intersection: &north_pole,
            },
            Test {
                rect: &rect,
                other: &Rect::from_degrees(-10., -1., 1., 20.),
                contains: false,
                intersects: true,
                union: &Rect::from_degrees(-10., 180., 90., 20.),
                intersection: &Rect::from_degrees(0., -1., 1., 0.),
            },
            Test {
                rect: &rect,
                other: &Rect::from_degrees(-10., -1., 0., 20.),
                contains: false,
                intersects: true,
                union: &Rect::from_degrees(-10., 180., 90., 20.),
                intersection: &Rect::from_degrees(0., -1., 0., 0.),
            },
            Test {
                rect: &rect,
                other: &Rect::from_degrees(-10., 0., 1., 20.),
                contains: false,
                intersects: true,
                union: &Rect::from_degrees(-10., 180., 90., 20.),
                intersection: &Rect::from_degrees(0., 0., 1., 0.),
            },
            Test {
                rect: &Rect::from_degrees(-15., -160., -15., -150.),
                other: &Rect::from_degrees(20., 145., 25., 155.),
                contains: false,
                intersects: false,
                union: &Rect::from_degrees(-15., 145., 25., -150.),
                intersection: &Rect::empty(),
            },
            Test {
                rect: &Rect::from_degrees(70., -10., 90., -140.),
                other: &Rect::from_degrees(60., 175., 80., 5.),
                contains: false,
                intersects: true,
                union: &Rect::from_degrees(60., -180., 90., 180.),
                intersection: &Rect::from_degrees(70., 175., 80., 5.),
            },
            // Check that the intersection of two rectangles that overlap in latitude
            // but not longitude is valid, and vice versa.
            Test {
                rect: &Rect::from_degrees(12., 30., 60., 60.),
                other: &Rect::from_degrees(0., 0., 30., 18.),
                contains: false,
                intersects: false,
                union: &Rect::from_degrees(0., 0., 60., 60.),
                intersection: &Rect::empty(),
            },
            Test {
                rect: &Rect::from_degrees(0., 0., 18., 42.),
                other: &Rect::from_degrees(30., 12., 42., 60.),
                contains: false,
                intersects: false,
                union: &Rect::from_degrees(0., 0., 42., 60.),
                intersection: &Rect::empty(),
            },
        ];
        for test in &tests {
            assert_eq!(test.rect.contains(test.other), test.contains);

            assert_eq!(test.rect.intersects(test.other), test.intersects);

            assert_eq!(
                test.rect.approx_eq(&test.rect.union(test.other)),
                test.rect.contains(test.other)
            );

            assert_ne!(
                test.rect.intersection(test.other).is_empty(),
                test.rect.intersects(test.other)
            );

            assert!(test.union.approx_eq(&test.rect.union(test.other)));

            assert!(test
                .intersection
                .approx_eq(&test.rect.intersection(test.other)));
        }
    }

    #[test]
    fn test_rect_cell_ops() {
        let cell0 = Cell::from(Point(Vector {
            x: 1. + 1e-12,
            y: 1.,
            z: 1.,
        }));
        let v0 = LatLng::from(cell0.vertex(0));

        let cell202 = Cell::from(CellID::from_face_pos_level(2, 0, 2));
        let bound202 = cell202.rect_bound();

        struct Test<'a> {
            r: &'a Rect,
            c: &'a Cell,
            contains: bool,
            intersects: bool,
        }
        let tests: [Test; 16] = [
            // Special cases
            Test {
                r: &Rect::empty(),
                c: &Cell::from(CellID::from_face_pos_level(3, 0, 0)),
                contains: false,
                intersects: false,
            },
            Test {
                r: &Rect::full(),
                c: &Cell::from(CellID::from_face_pos_level(2, 0, 0)),
                contains: true,
                intersects: true,
            },
            Test {
                r: &Rect::full(),
                c: &Cell::from(CellID::from_face_pos_level(5, 0, 25)),
                contains: true,
                intersects: true,
            },
            // This rectangle includes the first quadrant of face 0.  It's expanded
            // slightly because cell bounding rectangles are slightly conservative.
            Test {
                r: &Rect::from_degrees(-45.1, -45.1, 0.1, 0.1),
                c: &Cell::from(CellID::from_face_pos_level(0, 0, 0)),
                contains: false,
                intersects: true,
            },
            Test {
                r: &Rect::from_degrees(-45.1, -45.1, 0.1, 0.1),
                c: &Cell::from(CellID::from_face_pos_level(0, 0, 1)),
                contains: true,
                intersects: true,
            },
            Test {
                r: &Rect::from_degrees(-45.1, -45.1, 0.1, 0.1),
                c: &Cell::from(CellID::from_face_pos_level(1, 0, 1)),
                contains: false,
                intersects: false,
            },
            // This rectangle intersects the first quadrant of face 0.
            Test {
                r: &Rect::from_degrees(-10., -45., 10., 0.),
                c: &Cell::from(CellID::from_face_pos_level(0, 0, 0)),
                contains: false,
                intersects: true,
            },
            Test {
                r: &Rect::from_degrees(-10., -45., 10., 0.),
                c: &Cell::from(CellID::from_face_pos_level(0, 0, 1)),
                contains: false,
                intersects: true,
            },
            Test {
                r: &Rect::from_degrees(-10., -45., 10., 0.),
                c: &Cell::from(CellID::from_face_pos_level(1, 0, 1)),
                contains: false,
                intersects: false,
            },
            // Rectangle consisting of a single point.
            Test {
                r: &Rect::from_degrees(4., 4., 4., 4.),
                c: &Cell::from(CellID::from_face(0)),
                contains: false,
                intersects: true,
            },
            // Rectangles that intersect the bounding rectangle of a face
            // but not the face itself.
            Test {
                r: &Rect::from_degrees(41., -87., 42., -79.),
                c: &Cell::from(CellID::from_face(2)),
                contains: false,
                intersects: false,
            },
            Test {
                r: &Rect::from_degrees(-41., 160., -40., -160.),
                c: &Cell::from(CellID::from_face(5)),
                contains: false,
                intersects: false,
            },
            Test {
                // This is the leaf cell at the top right hand corner of face 0.
                // It has two angles of 60 degrees and two of 120 degrees.
                r: &Rect::from_degrees(
                    v0.lat.deg() - 1e-8,
                    v0.lng.deg() - 1e-8,
                    v0.lat.deg() - 2e-10,
                    v0.lng.deg() + 1e-10,
                ),
                c: &cell0,
                contains: false,
                intersects: false,
            },
            Test {
                // Rectangles that intersect a face but where no vertex of one region
                // is contained by the other region.  The first one passes through
                // a corner of one of the face cells.
                r: &Rect::from_degrees(-37., -70., -36., -20.),
                c: &Cell::from(CellID::from_face(5)),
                contains: false,
                intersects: true,
            },
            Test {
                // These two intersect like a diamond and a square.
                r: &Rect::from_degrees(
                    bound202.lo().lat.deg() + 3.,
                    bound202.lo().lng.deg() + 3.,
                    bound202.hi().lat.deg() - 3.,
                    bound202.hi().lng.deg() - 3.,
                ),
                c: &cell202,
                contains: false,
                intersects: true,
            },
            Test {
                // from a bug report
                r: &Rect::from_degrees(34.2572864, 135.2673642, 34.2707907, 135.2995742),
                c: &Cell::from(CellID(0x6007500000000000)),
                contains: false,
                intersects: true,
            },
        ];

        for test in &tests {
            assert_eq!(test.r.contains_cell(test.c), test.contains);
            assert_eq!(test.r.intersects_cell(test.c), test.intersects);
        }
    }

    #[test]
    fn test_rect_contains_point() {
        let r1 = Rect::from_degrees(0., -180., 90., 0.);

        let tests = [
            (
                &r1,
                &Point(Vector {
                    x: 0.5,
                    y: -0.3,
                    z: 0.1,
                }),
                true,
            ),
            (
                &r1,
                &Point(Vector {
                    x: 0.5,
                    y: 0.2,
                    z: 0.1,
                }),
                false,
            ),
        ];
        for &(r, p, want) in &tests {
            assert_eq!(r.contains_point(p), want);
        }
    }

    #[test]
    fn test_rect_intersects_lat_edge() {
        let tests = [
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.,
                }),
                Angle::from(Deg(41.)),
                Angle::from(Deg(-87.)).rad(),
                Angle::from(Deg(-79.)).rad(),
                false,
            ),
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.,
                }),
                Angle::from(Deg(42.)),
                Angle::from(Deg(-87.)).rad(),
                Angle::from(Deg(-79.)).rad(),
                false,
            ),
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: -1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                }),
                Angle::from(Deg(-3.)),
                Angle::from(Deg(-1.)).rad(),
                Angle::from(Deg(23.)).rad(),
                false,
            ),
            (
                Point(Vector {
                    x: 1.,
                    y: 0.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 0.,
                }),
                Angle::from(Deg(-28.)),
                Angle::from(Deg(69.)).rad(),
                Angle::from(Deg(115.)).rad(),
                false,
            ),
            (
                Point(Vector {
                    x: 0.,
                    y: 1.,
                    z: 0.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: -1.,
                }),
                Angle::from(Deg(44.)),
                Angle::from(Deg(60.)).rad(),
                Angle::from(Deg(177.)).rad(),
                false,
            ),
            (
                Point(Vector {
                    x: 0.,
                    y: 1.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 0.,
                    y: 1.,
                    z: -1.,
                }),
                Angle::from(Deg(-25.)),
                Angle::from(Deg(-74.)).rad(),
                Angle::from(Deg(-165.)).rad(),
                true,
            ),
            (
                Point(Vector {
                    x: 1.,
                    y: 0.,
                    z: 0.,
                }),
                Point(Vector {
                    x: 0.,
                    y: 0.,
                    z: -1.,
                }),
                Angle::from(Deg(-4.)),
                Angle::from(Deg(-152.)).rad(),
                Angle::from(Deg(171.)).rad(),
                true,
            ),
            // from a bug report
            (
                Point(Vector {
                    x: -0.589375791872893683986945,
                    y: 0.583248451588733285433364,
                    z: 0.558978908075738245564423,
                }),
                Point(Vector {
                    x: -0.587388131301997518107783,
                    y: 0.581281455376392863776402,
                    z: 0.563104832905072516524569,
                }),
                Angle::from(Deg(34.2572864)),
                2.3608609,
                2.3614230,
                true,
            ),
        ];

        for (a, b, lat, lng_lo, lng_hi, want) in &tests {
            assert_eq!(
                intersects_lat_edge(a, b, *lat, Interval::new(*lng_lo, *lng_hi)),
                *want
            );
        }
    }

    #[test]
    fn test_rect_intersects_lng_edge() {
        let tests = [
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.,
                }),
                Angle::from(Deg(41.)).rad(),
                Angle::from(Deg(42.)).rad(),
                Angle::from(Deg(-79.)),
                false,
            ),
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.,
                }),
                Angle::from(Deg(41.)).rad(),
                Angle::from(Deg(42.)).rad(),
                Angle::from(Deg(-87.)),
                false,
            ),
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.,
                }),
                Angle::from(Deg(42.)).rad(),
                Angle::from(Deg(41.)).rad(),
                Angle::from(Deg(79.)),
                false,
            ),
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: 1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.,
                }),
                Angle::from(Deg(41.)).rad(),
                Angle::from(Deg(42.)).rad(),
                Angle::from(Deg(87.)),
                false,
            ),
            (
                Point(Vector {
                    x: 0.,
                    y: -1.,
                    z: -1.,
                }),
                Point(Vector {
                    x: -1.,
                    y: 0.,
                    z: -1.,
                }),
                Angle::from(Deg(-87.)).rad(),
                Angle::from(Deg(13.)).rad(),
                Angle::from(Deg(-143.)),
                true,
            ),
            (
                Point(Vector {
                    x: 1.,
                    y: 1.,
                    z: -1.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.,
                }),
                Angle::from(Deg(-64.)).rad(),
                Angle::from(Deg(13.)).rad(),
                Angle::from(Deg(40.)),
                true,
            ),
            (
                Point(Vector {
                    x: 1.,
                    y: 1.,
                    z: 0.,
                }),
                Point(Vector {
                    x: -1.,
                    y: 0.,
                    z: -1.,
                }),
                Angle::from(Deg(-64.)).rad(),
                Angle::from(Deg(56.)).rad(),
                Angle::from(Deg(151.)),
                true,
            ),
            (
                Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: 0.,
                }),
                Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: -1.,
                }),
                Angle::from(Deg(-50.)).rad(),
                Angle::from(Deg(18.)).rad(),
                Angle::from(Deg(-84.)),
                true,
            ),
        ];

        for (a, b, lat_lo, lat_hi, lng, want) in &tests {
            assert_eq!(
                intersects_lng_edge(
                    a,
                    b,
                    r1::interval::Interval {
                        lo: *lat_lo,
                        hi: *lat_hi
                    },
                    *lng
                ),
                *want
            );
        }
    }

    // interval_distance returns the minimum distance (in radians) from X to the latitude
    // line segment defined by the given latitude and longitude interval.
    fn interval_distance(x: &LatLng, lat: Angle, iv: Interval) -> Angle {
        // Is x inside the longitude interval?
        if iv.contains(x.lng.rad()) {
            return Rad((x.lat.rad() - lat.rad()).abs()).into(); //FIXME
        }

        return x
            .distance(&LatLng {
                lat,
                lng: Rad(iv.lo).into(),
            })
            .min(x.distance(&LatLng {
                lat,
                lng: Rad(iv.hi).into(),
            }));
    }

    // Returns the minimum distance from X to the latitude line segment defined by
    // the given latitude and longitude interval.
    fn brute_force_rect_latlng_distance(r: &Rect, ll: &LatLng) -> Angle {
        let pt = &Point::from(ll);
        if r.contains_point(pt) {
            return Rad(0.).into();
        }

        let lo_lat = interval_distance(ll, Rad(r.lat.lo).into(), r.lng);
        let hi_lat = interval_distance(ll, Rad(r.lat.hi).into(), r.lng);
        let lo_lng = distance_from_segment(
            &Point::from(ll),
            &Point::from(LatLng {
                lat: Rad(r.lat.lo).into(),
                lng: Rad(r.lng.lo).into(),
            }),
            &Point::from(LatLng {
                lat: Rad(r.lat.hi).into(),
                lng: Rad(r.lng.lo).into(),
            }),
        );
        let hi_lng = distance_from_segment(
            &Point::from(ll),
            &Point::from(LatLng {
                lat: Rad(r.lat.lo).into(),
                lng: Rad(r.lng.hi).into(),
            }),
            &Point::from(LatLng {
                lat: Rad(r.lat.hi).into(),
                lng: Rad(r.lng.hi).into(),
            }),
        );

        return lo_lat.min(hi_lat).min(lo_lng).min(hi_lng);
    }

    #[test]
    fn test_distance_rect_from_latlng() {
        // Rect that spans 180.
        let a = &Rect::from(LatLng::from_degrees(-1., -1.)).add(&LatLng::from_degrees(2., 1.));
        // Rect near north pole.
        let b = &Rect::from(LatLng::from_degrees(86., 0.)).add(&LatLng::from_degrees(88., 2.));
        // Rect that touches north pole.
        let c = &Rect::from(LatLng::from_degrees(88., 0.)).add(&LatLng::from_degrees(90., 2.));

        let tests = [
            (a, -2., -1.),
            (a, 1., 2.),
            (b, 87., 3.),
            (b, 87., -1.),
            (b, 89., 1.),
            (b, 89., 181.),
            (b, 85., 1.),
            (b, 85., 181.),
            (b, 90., 0.),
            (c, 89., 3.),
            (c, 89., 90.),
            (c, 89., 181.),
        ];

        for &(r, lat, lng) in &tests {
            let ll = &LatLng::from_degrees(lat, lng);
            let got = r.distance_to_latlng(ll);
            let want = brute_force_rect_latlng_distance(r, ll);
            assert!(f64_near(got.rad(), want.rad(), 1e-10));
        }
    }

    #[test]
    fn test_distance_rect_from_latlng_random_pairs() {
        let mut rng = random::rng();
        for _ in 0..10000 {
            let r = random::rect(&mut rng);
            let ll = random::latlng(&mut rng);
            let got = r.distance_to_latlng(&ll);
            let want = brute_force_rect_latlng_distance(&r, &ll);
            assert!(f64_near(got.rad(), want.rad(), 1e-10));
        }
    }

    // This function assumes that DirectedHausdorffDistance() always
    // returns a distance from some point in a to b. So the function
    // mainly tests whether the returned distance is large enough,
    // and only does a weak test on whether it is small enough.
    fn verify_directed_hausdorff_distance(a: &Rect, b: &Rect) {
        let resolution = 0.1;
        let mut max_distance = Angle::default();
        let sample_size_on_lat = (a.lat.len() / resolution).round() as i32 + 1;
        let sample_size_on_lng = (a.lng.len() / resolution).round() as i32 + 1;
        let delta_on_lat = a.lat.len() / sample_size_on_lat as f64;
        let delta_on_lng = a.lng.len() / sample_size_on_lng as f64;
        let mut ll = LatLng::new(Angle::default(), Angle::from(Rad(a.lng.lo)));
        for _ in 0..=sample_size_on_lng {
            ll.lat = Angle::from(Rad(a.lat.lo));
            for _ in 0..=sample_size_on_lat {
                let d = b.distance_to_latlng(&ll.normalized());
                max_distance = max_distance.max(d);
                ll.lat = ll.lat + delta_on_lat;
            }
            ll.lng = ll.lng + delta_on_lng;
        }

        let got = a.directed_hausdorff_distance(b);
        assert!(
            max_distance.rad() <= got.rad() + 1e-10,
            "hausdorff({:?}, {:?}) = {:?} < {:?} - ε, but shouldn’t",
            a,
            b,
            got,
            max_distance
        );
        assert!(
            max_distance.rad() >= got.rad() - resolution,
            "hausdorff({:?}, {:?}) = {:?} > {:?} + resolution, but shouldn’t",
            a,
            b,
            got,
            max_distance
        );
    }

    #[test]
    fn test_directed_hausdorff_distance_random_pairs() {
        let mut rng = random::rng();
        for _ in 0..1000 {
            let a = random::rect(&mut rng);
            let b = random::rect(&mut rng);

            // a and b are *minimum* bounding rectangles of two random points,
            // in particular, their Voronoi diagrams are always of the same
            // topology. We take the "complements" of a and b for more
            // thorough testing.
            let a2 = Rect {
                lat: a.lat,
                lng: a.lng.complement(),
            };
            let b2 = Rect {
                lat: b.lat,
                lng: b.lng.complement(),
            };

            // Note that "a" and "b" come from the same distribution, so
            // there is no need to test pairs such as (b, a), (b, a2), etc.
            verify_directed_hausdorff_distance(&a, &b);
            verify_directed_hausdorff_distance(&a2, &b);
            verify_directed_hausdorff_distance(&a, &b2);
            verify_directed_hausdorff_distance(&a2, &b2);
        }
    }

    #[test]
    fn test_directed_hausdorff_distance_contained() {
        // Caller rect is contained in callee rect. Should return 0.
        let a = Rect::from_degrees(-10., 20., -5., 90.);
        let tests = vec![
            Rect::from_degrees(-10., 20., -5., 90.),
            Rect::from_degrees(-10., 19., -5., 91.),
            Rect::from_degrees(-11., 20., -4., 90.),
            Rect::from_degrees(-11., 19., -4., 91.),
        ];
        for test in tests {
            let got = a.directed_hausdorff_distance(&test);
            assert_f64_eq!(got.rad(), 0.0);
        }
    }

    #[test]
    fn test_directed_hausdorff_distance_point_to_rect() {
        // The Hausdorff distance from a point to a rect should be
        // the same as its distance to the rect.
        let a1 = LatLng::from_degrees(5., 8.);
        let a2 = LatLng::from_degrees(90., 10.); // North pole.
        let tests = vec![
            (a1, Rect::from_degrees(-85., -50., -80., 10.)),
            (a2, Rect::from_degrees(-85., -50., -80., 10.)),
            (a1, Rect::from_degrees(4., -10., 80., 10.)),
            (a2, Rect::from_degrees(4., -10., 80., 10.)),
            (a1, Rect::from_degrees(70., 170., 80., -170.)),
            (a2, Rect::from_degrees(70., 170., 80., -170.)),
        ];
        for test in tests {
            let a = Rect::from(test.0);
            let got = a.directed_hausdorff_distance(&test.1).rad();
            let want = test.1.distance_to_latlng(&test.0).rad();
            assert_f64_eq!(got, want);
        }
    }

    #[test]
    fn test_directed_hausdorff_distance_rect_to_point() {
        let a = Rect::from_degrees(1., -8., 10., 20.);
        let tests = vec![(5, 8), (-6, -100), (-90, -20), (90, 0)];
        for test in tests {
            let ll = LatLng::from_degrees(test.0 as f64, test.1 as f64);
            verify_directed_hausdorff_distance(&a, &Rect::from(ll));
        }
    }

    #[test]
    fn test_directed_hausdorff_distance_rect_to_rect_near_pole() {
        let a = Rect::from_degrees(-87., 0., -85., 3.);
        let tests = vec![
            Rect::from_degrees(-89., 1., -88., 2.),
            Rect::from_degrees(-84., 1., -83., 2.),
            Rect::from_degrees(-88., 90., -86., 91.),
            Rect::from_degrees(-84., -91., -83., -90.),
            Rect::from_degrees(-90., 181., -89., 182.),
            Rect::from_degrees(-84., 181., -83., 182.),
        ];
        for test in tests {
            verify_directed_hausdorff_distance(&a, &test);
        }
    }

    #[test]
    fn test_directed_hausdorff_distance_rect_to_rect_degenerate_cases() {
        // Rectangles that contain poles.
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(0., 10., 90., 20.),
            &Rect::from_degrees(-4., -10., 4., 0.),
        );
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(-4., -10., 4., 0.),
            &Rect::from_degrees(0., 10., 90., 20.),
        );

        // Two rectangles share same or complement longitudinal intervals.
        let a = Rect::from_degrees(-50., -10., 50., 10.);
        let b = Rect::from_degrees(30., -10., 60., 10.);
        verify_directed_hausdorff_distance(&a, &b);

        let c = Rect {
            lat: a.lat,
            lng: a.lng.complement(),
        };
        verify_directed_hausdorff_distance(&c, &b);

        // Rectangle a touches b_opposite_lng.
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(10., 170., 30., 180.),
            &Rect::from_degrees(-50., -10., 50., 10.),
        );
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(10., -180., 30., -170.),
            &Rect::from_degrees(-50., -10., 50., 10.),
        );

        // Rectangle b's Voronoi diagram is degenerate (lng interval spans
        // 180 degrees), and a touches the degenerate Voronoi vertex.
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(-30., 170., 30., 180.),
            &Rect::from_degrees(-10., -90., 10., 90.),
        );
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(-30., -180., 30., -170.),
            &Rect::from_degrees(-10., -90., 10., 90.),
        );

        // Rectangle a touches a Voronoi vertex of rectangle b.
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(-20., 105., 20., 110.),
            &Rect::from_degrees(-30., 5., 30., 15.),
        );
        verify_directed_hausdorff_distance(
            &Rect::from_degrees(-20., 95., 20., 105.),
            &Rect::from_degrees(-30., 5., 30., 15.),
        );
    }

    #[test]
    fn test_rect_approx_equal() {
        // s1.Interval and r1.Interval have additional testing.

        let e: f64 = 1e-15;
        let tests = [
            (Rect::empty(), Rect::from_degrees(1., 5., 1., 5.), true),
            (Rect::from_degrees(1., 5., 1., 5.), Rect::empty(), true),
            (
                Rect::from_degrees(1., 5., 1., 5.),
                Rect::from_degrees(2., 7., 2., 7.),
                false,
            ),
            (
                Rect::from_degrees(1., 5., 1., 5.),
                Rect::from_degrees(1. + e, 5. + e, 1. + e, 5. + e),
                true,
            ),
        ];

        for (a, b, want) in &tests {
            assert_eq!(a.approx_eq(b), *want);
        }
    }

    #[test]
    fn test_rect_centroid_empty_full() {
        assert_eq!(Rect::empty().centroid(), Point::default());
        assert_f64_eq!(Rect::full().centroid().norm(), EPSILON);
    }

    /// Recursively verify that when a rectangle is split into two pieces,
    /// the centroids of the children sum to give the centroid of the parent.
    fn test_rect_centroid_splitting(r: Rect, splits_left: i32) {
        let mut rng = random::rng();
        let child0: Rect;
        let child1: Rect;
        if random::one_in(&mut rng, 2) {
            let lat_lo = r.lat.lo.min(r.lat.hi);
            let lat_hi = r.lat.lo.max(r.lat.hi);
            let lat = rng.gen_range(lat_lo..lat_hi);
            child0 = Rect {
                lat: r1::interval::Interval::new(r.lat.lo, lat),
                lng: r.lng,
            };
            child1 = Rect {
                lat: r1::interval::Interval::new(lat, r.lat.hi),
                lng: r.lng,
            };
        } else {
            let lng_lo = r.lng.lo.min(r.lng.hi);
            let lng_hi = r.lng.lo.max(r.lng.hi);
            let lng = rng.gen_range(lng_lo..lng_hi);
            child0 = Rect {
                lat: r.lat,
                lng: Interval {
                    lo: r.lng.lo,
                    hi: lng,
                },
            };
            child1 = Rect {
                lat: r.lat,
                lng: Interval {
                    lo: lng,
                    hi: r.lng.hi,
                },
            };
        }
        let got = (r.centroid() - child0.centroid() - child1.centroid()).norm();
        assert!(
            got <= EPSILON,
            "want ~0, got {:?}, r={:?}, child0={:?}, child1={:?}",
            got,
            r,
            child0,
            child1
        );
        if splits_left > 0 {
            test_rect_centroid_splitting(child0, splits_left - 1);
            test_rect_centroid_splitting(child1, splits_left - 1);
        }
    }

    #[test]
    fn test_rect_centroid_full_range() {
        let mut rng = random::rng();

        // Rectangles that cover the full longitude range.
        for _ in 0..100 {
            let lat1 = rng.gen_range(-FRAC_PI_2..FRAC_PI_2);
            let lat2 = rng.gen_range(-FRAC_PI_2..FRAC_PI_2);
            let r = Rect {
                lat: r1::interval::Interval::new(lat1, lat2),
                lng: VALID_RECT_LNG_RANGE,
            };
            let centroid = r.centroid();
            let want = 0.5 * (lat1.sin() + lat2.sin()) * r.area();
            assert_f64_eq!(centroid.0.z, want);
        }

        // Rectangles that cover the full latitude range.
        for _ in 0..100 {
            let lng1 = rng.gen_range(-PI..PI);
            let lng2 = rng.gen_range(-PI..PI);
            let r = Rect {
                lat: VALID_RECT_LAT_RANGE,
                lng: Interval::new(lng1, lng2),
            };
            let centroid = r.centroid();
            assert!(centroid.0.z.abs() <= EPSILON);
            assert_f64_eq!(r.lng.center(), LatLng::from(centroid).lng.rad());
            let alpha = 0.5 * r.lng.len();
            let p = crate::r2::point::Point::new(centroid.0.x, centroid.0.y);
            assert_f64_eq!(0.25 * PI * alpha.sin() / alpha * r.area(), p.norm());
        }

        // Finally, verify that when a rectangle is recursively split
        // into pieces, the centroids of the pieces add to give the
        // centroid of their parent.  To make the code simpler we avoid
        // rectangles that cross the 180 degree line of longitude.
        test_rect_centroid_splitting(
            Rect {
                lat: VALID_RECT_LAT_RANGE,
                lng: VALID_RECT_LNG_RANGE,
            },
            10,
        );
    }
}
