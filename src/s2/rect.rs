use std;
use std::f64::consts::PI;

use crate::consts::*;
use crate::r1;
use crate::s1::*;
use crate::s2::edgeutil;
use crate::s2::latlng::LatLng;

#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Rect {
    pub lat: r1::interval::Interval,
    pub lng: Interval,
}

impl std::fmt::Debug for Rect {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[lo{:?}, hi{:?}]", self.lo(), self.hi())
    }
}

const VALID_RECT_LAT_RANGE: r1::interval::Interval = r1::interval::Interval {
    lo: -PI / 2.,
    hi: PI / 2.,
};
const VALID_RECT_LNG_RANGE: Interval = interval::FULL;

impl Rect {
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
            lat: r1::interval::Interval {
                lo: Angle::from(Deg(lat_lo)).rad(),
                hi: Angle::from(Deg(lat_hi)).rad(),
            },
            lng: Interval::new(
                Angle::from(Deg(lng_lo)).rad(),
                Angle::from(Deg(lng_hi)).rad(),
            ),
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

    pub fn vertex(&self, i: u8) -> LatLng {
        let (lat, lng) = match i {
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

    // extra functions
    pub fn approx_eq(&self, other: &Self) -> bool {
        self.lat.approx_eq(&other.lat) && self.lng.approx_eq(&other.lng)
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

    fn intersects_cell(&self, cell: &Cell) -> bool {
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
        let mut vertices = Vec::new();
        let mut latlngs = Vec::new();

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
            let edge_lng =
                interval::Interval::new(latlngs[i].lng.rad(), latlngs[(i + 1) & 3].lng.rad());
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
fn intersects_lat_edge(a: &Point, b: &Point, lat: Angle, lng: interval::Interval) -> bool {
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
    let ab_theta = interval::Interval::from_point_pair(
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
}

/*

// IntersectsCell reports whether this rectangle intersects the given cell. This is an
// exact test and may be fairly expensive.
func (r Rect) IntersectsCell(c Cell) bool {
    // First we eliminate the cases where one region completely contains the
    // other. Once these are disposed of, then the regions will intersect
    // if and only if their boundaries intersect.
    if r.IsEmpty() {
        return false
    }
    if r.ContainsPoint(Point{c.id.rawPoint()}) {
        return true
    }
    if c.ContainsPoint(PointFromLatLng(r.Center())) {
        return true
    }

    // Quick rejection test (not required for correctness).
    if !r.Intersects(c.RectBound()) {
        return false
    }

    // Precompute the cell vertices as points and latitude-longitudes. We also
    // check whether the Cell contains any corner of the rectangle, or
    // vice-versa, since the edge-crossing tests only check the edge interiors.
    vertices := [4]Point{}
    latlngs := [4]LatLng{}

    for i := range vertices {
        vertices[i] = c.Vertex(i)
        latlngs[i] = LatLngFromPoint(vertices[i])
        if r.ContainsLatLng(latlngs[i]) {
            return true
        }
        if c.ContainsPoint(PointFromLatLng(r.Vertex(i))) {
            return true
        }
    }

    // Now check whether the boundaries intersect. Unfortunately, a
    // latitude-longitude rectangle does not have straight edges: two edges
    // are curved, and at least one of them is concave.
    for i := range vertices {
        edgeLng := s1.IntervalFromEndpoints(latlngs[i].Lng.Radians(), latlngs[(i+1)&3].Lng.Radians())
        if !r.Lng.Intersects(edgeLng) {
            continue
        }

        a := vertices[i]
        b := vertices[(i+1)&3]
        if edgeLng.Contains(r.Lng.Lo) && intersectsLngEdge(a, b, r.Lat, s1.Angle(r.Lng.Lo)) {
            return true
        }
        if edgeLng.Contains(r.Lng.Hi) && intersectsLngEdge(a, b, r.Lat, s1.Angle(r.Lng.Hi)) {
            return true
        }
        if intersectsLatEdge(a, b, s1.Angle(r.Lat.Lo), r.Lng) {
            return true
        }
        if intersectsLatEdge(a, b, s1.Angle(r.Lat.Hi), r.Lng) {
            return true
        }
    }
    return false
}

// BUG: The major differences from the C++ version are:
//   - GetCentroid, Get*Distance, Vertex, InteriorContains(LatLng|Rect|Point)
*/

#[cfg(test)]
#[allow(non_upper_case_globals)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    use crate::cellid::CellID;
    use crate::predicates::sign;
    use crate::r1;
    use crate::r3::vector::Vector;
    use std::ops::Add;

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

    // TODO: test_rect_string

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
    fn test_rect_vertex() {
        let r1 = &Rect {
            lat: r1::interval::Interval {
                lo: 0.0,
                hi: PI / 2.0,
            },
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
        };
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
                lat: lat,
                lng: Rad(iv.lo).into(),
            })
            .min(x.distance(&LatLng {
                lat: lat,
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
    /*
        func TestDistanceRectFromLatLngRandomPairs(t *testing.T) {
            latlng := func() LatLng { return LatLngFromPoint(randomPoint()) }

            for i := 0; i < 10000; i++ {
                r := RectFromLatLng(latlng()).AddPoint(latlng())
                ll := latlng()
                got := r.DistanceToLatLng(ll)
                want := bruteForceRectLatLngDistance(r, ll)
                if !float64Near(float64(got), float64(want), 1e-10) {
                    t.Errorf("dist from %v to %v = %v, want %v", r, ll, got, want)
                }
            }
        }

        // This function assumes that DirectedHausdorffDistance() always returns
        // a distance from some point in a to b. So the function mainly tests whether
        // the returned distance is large enough, and only does a weak test on whether
        // it is small enough.
        func verifyDirectedHausdorffDistance(t *testing.T, a, b Rect) {
            t.Helper()

            const resolution = 0.1

            // Record the max sample distance as well as the sample point realizing the
            // max for easier debugging.
            var maxDistance s1.Angle

            sampleSizeOnLat := int(a.Lat.Length()/resolution) + 1
            sampleSizeOnLng := int(a.Lng.Length()/resolution) + 1

            deltaOnLat := s1.Angle(a.Lat.Length()) / s1.Angle(sampleSizeOnLat)
            deltaOnLng := s1.Angle(a.Lng.Length()) / s1.Angle(sampleSizeOnLng)

            ll := LatLng{Lng: s1.Angle(a.Lng.Lo)}
            for i := 0; i <= sampleSizeOnLng; i++ {
                ll.Lat = s1.Angle(a.Lat.Lo)

                for j := 0; j <= sampleSizeOnLat; j++ {
                    d := b.DistanceToLatLng(ll.Normalized())
                    maxDistance = maxAngle(maxDistance, d)
                    ll.Lat += deltaOnLat
                }
                ll.Lng += deltaOnLng
            }

            got := a.DirectedHausdorffDistance(b)

            if got < maxDistance-1e-10 {
                t.Errorf("hausdorff(%v, %v) = %v < %v-eps, but shouldn't", a, b, got, maxDistance)
            } else if got > maxDistance+resolution {
                t.Errorf("DirectedHausdorffDistance(%v, %v) = %v > %v+resolution, but shouldn't", a, b, got, maxDistance)
            }
        }

        func TestRectDirectedHausdorffDistanceRandomPairs(t *testing.T) {
            // Test random pairs.
            rnd := func() LatLng { return LatLngFromPoint(randomPoint()) }
            for i := 0; i < 1000; i++ {
                a := RectFromLatLng(rnd()).AddPoint(rnd())
                b := RectFromLatLng(rnd()).AddPoint(rnd())
                // a and b are *minimum* bounding rectangles of two random points, in
                // particular, their Voronoi diagrams are always of the same topology. We
                // take the "complements" of a and b for more thorough testing.
                a2 := Rect{Lat: a.Lat, Lng: a.Lng.Complement()}
                b2 := Rect{Lat: b.Lat, Lng: b.Lng.Complement()}

                // Note that "a" and "b" come from the same distribution, so there is no
                // need to test pairs such as (b, a), (b, a2), etc.
                verifyDirectedHausdorffDistance(t, a, b)
                verifyDirectedHausdorffDistance(t, a2, b)
                verifyDirectedHausdorffDistance(t, a, b2)
                verifyDirectedHausdorffDistance(t, a2, b2)
            }
        }

        func TestDirectedHausdorffDistanceContained(t *testing.T) {
            // Caller rect is contained in callee rect. Should return 0.
            a := rectFromDegrees(-10, 20, -5, 90)
            tests := []Rect{
                rectFromDegrees(-10, 20, -5, 90),
                rectFromDegrees(-10, 19, -5, 91),
                rectFromDegrees(-11, 20, -4, 90),
                rectFromDegrees(-11, 19, -4, 91),
            }
            for _, test := range tests {
                got, want := a.DirectedHausdorffDistance(test), s1.Angle(0)
                if got != want {
                    t.Errorf("%v.DirectedHausdorffDistance(%v) = %v, want %v", a, test, got, want)
                }
            }
        }

        func TestDirectHausdorffDistancePointToRect(t *testing.T) {
            // The Hausdorff distance from a point to a rect should be the same as its
            // distance to the rect.
            a1 := LatLngFromDegrees(5, 8)
            a2 := LatLngFromDegrees(90, 10) // North pole.

            tests := []struct {
                ll LatLng
                b  Rect
            }{
                {a1, rectFromDegrees(-85, -50, -80, 10)},
                {a2, rectFromDegrees(-85, -50, -80, 10)},
                {a1, rectFromDegrees(4, -10, 80, 10)},
                {a2, rectFromDegrees(4, -10, 80, 10)},
                {a1, rectFromDegrees(70, 170, 80, -170)},
                {a2, rectFromDegrees(70, 170, 80, -170)},
            }
            for _, test := range tests {
                a := RectFromLatLng(test.ll)
                got, want := a.DirectedHausdorffDistance(test.b), test.b.DistanceToLatLng(test.ll)

                if !float64Eq(float64(got), float64(want)) {
                    t.Errorf("hausdorff(%v, %v) = %v, want %v, as that's the closest dist", test.b, a, got, want)
                }
            }
        }

        func TestDirectedHausdorffDistanceRectToPoint(t *testing.T) {
            a := rectFromDegrees(1, -8, 10, 20)
            tests := []struct {
                lat, lng float64 // Degrees.
            }{{5, 8}, {-6, -100}, {-90, -20}, {90, 0}}
            for _, test := range tests {
                verifyDirectedHausdorffDistance(t, a, RectFromLatLng(LatLngFromDegrees(test.lat, test.lng)))
            }
        }

        func TestDirectedHausdorffDistanceRectToRectNearPole(t *testing.T) {
            // Tests near south pole.
            a := rectFromDegrees(-87, 0, -85, 3)
            tests := []Rect{
                rectFromDegrees(-89, 1, -88, 2),
                rectFromDegrees(-84, 1, -83, 2),
                rectFromDegrees(-88, 90, -86, 91),
                rectFromDegrees(-84, -91, -83, -90),
                rectFromDegrees(-90, 181, -89, 182),
                rectFromDegrees(-84, 181, -83, 182),
            }
            for _, test := range tests {
                verifyDirectedHausdorffDistance(t, a, test)
            }
        }

        func TestDirectedHausdorffDistanceRectToRectDegenerateCases(t *testing.T) {
            // Rectangles that contain poles.
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(0, 10, 90, 20), rectFromDegrees(-4, -10, 4, 0))
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(-4, -10, 4, 0), rectFromDegrees(0, 10, 90, 20))

            // Two rectangles share same or complement longitudinal intervals.
            a := rectFromDegrees(-50, -10, 50, 10)
            b := rectFromDegrees(30, -10, 60, 10)
            verifyDirectedHausdorffDistance(t, a, b)

            c := Rect{Lat: a.Lat, Lng: a.Lng.Complement()}
            verifyDirectedHausdorffDistance(t, c, b)

            // Rectangle a touches b_opposite_lng.
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(10, 170, 30, 180), rectFromDegrees(-50, -10, 50, 10))
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(10, -180, 30, -170), rectFromDegrees(-50, -10, 50, 10))

            // Rectangle b's Voronoi diagram is degenerate (lng interval spans 180
            // degrees), and a touches the degenerate Voronoi vertex.
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(-30, 170, 30, 180), rectFromDegrees(-10, -90, 10, 90))
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(-30, -180, 30, -170), rectFromDegrees(-10, -90, 10, 90))

            // Rectangle a touches a voronoi vertex of rectangle b.
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(-20, 105, 20, 110), rectFromDegrees(-30, 5, 30, 15))
            verifyDirectedHausdorffDistance(t,
                rectFromDegrees(-20, 95, 20, 105), rectFromDegrees(-30, 5, 30, 15))
        }
    */
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
    /*
    fn testRectCentroidSplitting() {
        // Recursively verify that when a rectangle is split into two pieces, the
        // centroids of the children sum to give the centroid of the parent.
        let xx = Rect::default();
        var child0, child1 Rect
        if oneIn(2) {
            lat := randomUniformFloat64(r.Lat.Lo, r.Lat.Hi)
            child0 = Rect{r1.Interval{r.Lat.Lo, lat}, r.Lng}
            child1 = Rect{r1.Interval{lat, r.Lat.Hi}, r.Lng}
        } else {
            lng := randomUniformFloat64(r.Lng.Lo, r.Lng.Hi)
            child0 = Rect{r.Lat, s1.Interval{r.Lng.Lo, lng}}
            child1 = Rect{r.Lat, s1.Interval{lng, r.Lng.Hi}}
        }

        if got, want := r.Centroid().Sub(child0.Centroid().Vector).Sub(child1.Centroid().Vector).Norm(), 1e-15; got > want {
            t.Errorf("%v.Centroid() - %v.Centroid() - %v.Centroid = %v, want ~0", r, child0, child1, got)
        }
        if leftSplits > 0 {
            testRectCentroidSplitting(t, child0, leftSplits-1)
            testRectCentroidSplitting(t, child1, leftSplits-1)
        }
    }*/
    /*
    func TestRectCentroidFullRange(t *testing.T) {
        // Rectangles that cover the full longitude range.
        for i := 0; i < 100; i++ {
            lat1 := randomUniformFloat64(-math.Pi/2, math.Pi/2)
            lat2 := randomUniformFloat64(-math.Pi/2, math.Pi/2)
            r := Rect{r1.Interval{lat1, lat2}, s1.FullInterval()}
            centroid := r.Centroid()
            if want := 0.5 * (math.Sin(lat1) + math.Sin(lat2)) * r.Area(); !float64Near(want, centroid.Z, epsilon) {
                t.Errorf("%v.Centroid().Z was %v, want %v", r, centroid.Z, want)
            }
            if got := (r2.Point{centroid.X, centroid.Y}.Norm()); got > epsilon {
                t.Errorf("%v.Centroid().Norm() was %v, want > %v ", r, got, epsilon)
            }
        }

        // Rectangles that cover the full latitude range.
        for i := 0; i < 100; i++ {
            lat1 := randomUniformFloat64(-math.Pi, math.Pi)
            lat2 := randomUniformFloat64(-math.Pi, math.Pi)
            r := Rect{r1.Interval{-math.Pi / 2, math.Pi / 2}, s1.Interval{lat1, lat2}}
            centroid := r.Centroid()

            if got, want := math.Abs(centroid.Z), epsilon; got > want {
                t.Errorf("math.Abs(%v.Centroid().Z) = %v, want <= %v", r, got, want)
            }

            if got, want := LatLngFromPoint(centroid).Lng.Radians(), r.Lng.Center(); !float64Near(got, want, epsilon) {
                t.Errorf("%v.Lng.Radians() = %v, want %v", centroid, got, want)
            }

            alpha := 0.5 * r.Lng.Length()
            if got, want := (r2.Point{centroid.X, centroid.Y}.Norm()), (0.25 * math.Pi * math.Sin(alpha) / alpha * r.Area()); !float64Near(got, want, epsilon) {
                t.Errorf("%v.Centroid().Norm() = %v, want ~%v", got, want, epsilon)
            }
        }

        // Finally, verify that when a rectangle is recursively split into pieces,
        // the centroids of the pieces add to give the centroid of their parent.
        // To make the code simpler we avoid rectangles that cross the 180 degree
        // line of longitude.
        testRectCentroidSplitting(t, Rect{r1.Interval{-math.Pi / 2, math.Pi / 2}, s1.Interval{-math.Pi, math.Pi}}, 10)
         */
}
