
use std;

use consts::*;
use r3::vector::Vector;
use s1;
use s1::chordangle::ChordAngle;
use s2::predicates::*;
use s2::region::Region;
use s2::latlng::LatLng;
use s2::rect::Rect;
use s2::cap::Cap;
use s2::cell::Cell;

/// Point represents a point on the unit sphere as a normalized 3D vector.
/// Fields should be treated as read-only. Use one of the factory methods for creation.
#[derive(Clone,PartialEq,Debug)]
pub struct Point(pub Vector);

impl std::ops::Add<Point> for Point {
    type Output = Point;
    fn add(self, other: Point) -> Self::Output {
        Point(self.0 + other.0)
    }
}

impl std::ops::Sub<Point> for Point {
    type Output = Point;
    fn sub(self, other: Point) -> Self::Output {
        &self - &other
    }
}
impl<'a, 'b> std::ops::Sub<&'b Point> for &'a Point {
    type Output = Point;
    fn sub(self, other: &'b Point) -> Self::Output {
        Point(&self.0 - &other.0)
    }
}

impl std::ops::Mul<Point> for Point {
    type Output = Point;
    fn mul(self, other: Point) -> Self::Output {
        Point(self.0 * other.0)
    }
}
impl std::ops::Mul<f64> for Point {
    type Output = Point;
    fn mul(self, m: f64) -> Self::Output {
        Point(self.0 * m)
    }
}

pub const ORIGIN: Point = Point(Vector {
                                    x: -0.0099994664350250197,
                                    y: 0.0025924542609324121,
                                    z: 0.99994664350250195,
                                });

impl Point {
    /// form_coords creates a new normalized point from coordinates.
    ///
    /// This always returns a valid point. If the given coordinates can not be normalized
    /// the origin point will be returned.
    ///
    /// This behavior is different from the C++ construction of a S2Point from coordinates
    /// (i.e. S2Point(x, y, z)) in that in C++ they do not Normalize.
    pub fn from_coords(x: f64, y: f64, z: f64) -> Self {
        if x == 0. && y == 0. && z == 0. {
            Point::origin()
        } else {
            Point(Vector { x: x, y: y, z: z }.normalize())
        }
    }

    /// origin returns a unique "origin" on the sphere for operations that need a fixed
    /// reference point. In particular, this is the "point at infinity" used for
    /// point-in-polygon testing (by counting the number of edge crossings).
    ///
    /// It should *not* be a point that is commonly used in edge tests in order
    /// to avoid triggering code to handle degenerate cases (this rules out the
    /// north and south poles). It should also not be on the boundary of any
    /// low-level S2Cell for the same reason.
    pub fn origin() -> Self {
        ORIGIN
    }

    /// cross returns a Point that is orthogonal to both p and op. This is similar to
    /// p.Cross(op) (the true cross product) except that it does a better job of
    /// ensuring orthogonality when the Point is nearly parallel to op, it returns
    /// a non-zero result even when p == op or p == -op and the result is a Point.
    ///
    /// It satisfies the following properties (f == cross):
    ///
    /// ```text
    /// (1) f(p, op) != 0 for all p, op
    /// (2) f(op,p) == -f(p,op) unless p == op or p == -op
    /// (3) f(-p,op) == -f(p,op) unless p == op or p == -op
    /// (4) f(p,-op) == -f(p,op) unless p == op or p == -op
    /// ```
    pub fn cross(&self, other: &Self) -> Self {
        // NOTE(dnadasi): In the C++ API the equivalent method here was known as "RobustCrossProd",
        // but point_cross more accurately describes how this method is used.
        let v = (&self.0 + &other.0).cross(&(&other.0 - &self.0));

        // Compare exactly to the 0 vector.
        if v.x == 0. && v.y == 0. && v.z == 0. {
            // The only result that makes sense mathematically is to return zero, but
            // we find it more convenient to return an arbitrary orthogonal vector.
            Point(self.0.ortho())
        } else {
            Point(v)
        }
    }

    /// distance returns the angle between two points.
    pub fn distance(&self, b: &Point) -> s1::angle::Angle {
        self.0.angle(&b.0)
    }

    /// approx_eq reports whether the two points are similar enough to be equal.
    pub fn approx_eq(&self, other: &Self) -> bool {
        self.0.angle(&other.0) <= s1::angle::Angle(EPSILON)
    }

    /// norm returns the point's norm.
    pub fn norm(&self) -> f64 {
        self.0.norm()
    }
}

/// ordered_ccw returns true if the edges OA, OB, and OC are encountered in that
/// order while sweeping CCW around the point O.
///
/// You can think of this as testing whether A <= B <= C with respect to the
/// CCW ordering around O that starts at A, or equivalently, whether B is
/// contained in the range of angles (inclusive) that starts at A and extends
/// CCW to C. Properties:
///
///  (1) If OrderedCCW(a,b,c,o) && OrderedCCW(b,a,c,o), then a == b
///  (2) If OrderedCCW(a,b,c,o) && OrderedCCW(a,c,b,o), then b == c
///  (3) If OrderedCCW(a,b,c,o) && OrderedCCW(c,b,a,o), then a == b == c
///  (4) If a == b or b == c, then OrderedCCW(a,b,c,o) is true
///  (5) Otherwise if a == c, then OrderedCCW(a,b,c,o) is false
pub fn ordered_ccw(a: &Point, b: &Point, c: &Point, o: &Point) -> bool {
    let mut sum = 0;
    if robust_sign(b, o, a) != Direction::Clockwise {
        sum += 1;
    }
    if robust_sign(c, o, b) != Direction::Clockwise {
        sum += 1;
    }
    if robust_sign(a, o, c) == Direction::CounterClockwise {
        sum += 1;
    }
    return sum >= 2;
}

/// point_area returns the area on the unit sphere for the triangle defined by the
/// given points.
///
/// This method is based on l'Huilier's theorem,
///
///   tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
///
/// where E is the spherical excess of the triangle (i.e. its area),
///       a, b, c are the side lengths, and
///       s is the semiperimeter (a + b + c) / 2.
///
/// The only significant source of error using l'Huilier's method is the
/// cancellation error of the terms (s-a), (s-b), (s-c). This leads to a
/// *relative* error of about 1e-16 * s / min(s-a, s-b, s-c). This compares
/// to a relative error of about 1e-15 / E using Girard's formula, where E is
/// the true area of the triangle. Girard's formula can be even worse than
/// this for very small triangles, e.g. a triangle with a true area of 1e-30
/// might evaluate to 1e-5.
///
/// So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
/// dmin = min(s-a, s-b, s-c). This basically includes all triangles
/// except for extremely long and skinny ones.
///
/// Since we don't know E, we would like a conservative upper bound on
/// the triangle area in terms of s and dmin. It's possible to show that
/// E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
/// Using this, it's easy to show that we should always use l'Huilier's
/// method if dmin >= k2 * s^5, where k2 is about 1e-2. Furthermore,
/// if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
/// k3 is about 0.1. Since the best case error using Girard's formula
/// is about 1e-15, this means that we shouldn't even consider it unless
/// s >= 3e-4 or so.
pub fn point_area(a: &Point, b: &Point, c: &Point) -> f64 {
    let sa = b.0.angle(&c.0).0;
    let sb = c.0.angle(&a.0).0;
    let sc = a.0.angle(&b.0).0;
    let s = 0.5 * (sa + sb + sc);
    if s >= 3e-4 {
        // Consider whether Girard's formula might be more accurate.
        let dmin = s - sa.max(sb.max(sc));
        if dmin < 1e-2 * s * s * s * s * s {
            // This triangle is skinny enough to use Girard's formula.
            let ab = a.cross(b);
            let bc = b.cross(c);
            let ac = a.cross(c);
            let area = (ab.0.angle(&ac.0).0 - ab.0.angle(&bc.0).0 + bc.0.angle(&ac.0).0).max(0.0);

            if dmin < s * 0.1 * area {
                return area;
            }
        }
    }

    // Use l'Huilier's formula.
    4. *
    ((0.5 * s).tan() * (0.5 * (s - sa)).tan() * (0.5 * (s - sb)).tan() * (0.5 * (s - sc)).tan())
        .max(0.)
        .sqrt()
        .atan()
}

/*
// TrueCentroid returns the true centroid of the spherical triangle ABC multiplied by the
// signed area of spherical triangle ABC. The result is not normalized.
// The reasons for multiplying by the signed area are (1) this is the quantity
// that needs to be summed to compute the centroid of a union or difference of triangles,
// and (2) it's actually easier to calculate this way. All points must have unit length.
//
// The true centroid (mass centroid) is defined as the surface integral
// over the spherical triangle of (x,y,z) divided by the triangle area.
// This is the point that the triangle would rotate around if it was
// spinning in empty space.
//
// The best centroid for most purposes is the true centroid. Unlike the
// planar and surface centroids, the true centroid behaves linearly as
// regions are added or subtracted. That is, if you split a triangle into
// pieces and compute the average of their centroids (weighted by triangle
// area), the result equals the centroid of the original triangle. This is
// not true of the other centroids.
func TrueCentroid(a, b, c Point) Point {
	ra := float64(1)
	if sa := float64(b.distance(c)); sa != 0 {
		ra = sa / math.Sin(sa)
	}
	rb := float64(1)
	if sb := float64(c.distance(a)); sb != 0 {
		rb = sb / math.Sin(sb)
	}
	rc := float64(1)
	if sc := float64(a.distance(b)); sc != 0 {
		rc = sc / math.Sin(sc)
	}

	// Now compute a point M such that:
	//
	//  [Ax Ay Az] [Mx]                       [ra]
	//  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
	//  [Cx Cy Cz] [Mz]                       [rc]
	//
	// To improve the numerical stability we subtract the first row (A) from the
	// other two rows; this reduces the cancellation error when A, B, and C are
	// very close together. Then we solve it using Cramer's rule.
	//
	// This code still isn't as numerically stable as it could be.
	// The biggest potential improvement is to compute B-A and C-A more
	// accurately so that (B-A)x(C-A) is always inside triangle ABC.
	x := r3.Vector{a.X, b.X - a.X, c.X - a.X}
	y := r3.Vector{a.Y, b.Y - a.Y, c.Y - a.Y}
	z := r3.Vector{a.Z, b.Z - a.Z, c.Z - a.Z}
	r := r3.Vector{ra, rb - ra, rc - ra}

	return Point{r3.Vector{y.Cross(z).Dot(r), z.Cross(x).Dot(r), x.Cross(y).Dot(r)}.Mul(0.5)}
}
*/

/// planar_centroid returns the centroid of the planar triangle ABC, which is not normalized.
/// It can be normalized to unit length to obtain the "surface centroid" of the corresponding
/// spherical triangle, i.e. the intersection of the three medians. However,
/// note that for large spherical triangles the surface centroid may be
/// nowhere near the intuitive "center" (see example in TrueCentroid comments).
///
/// Note that the surface centroid may be nowhere near the intuitive
/// "center" of a spherical triangle. For example, consider the triangle
/// with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
/// The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
/// within a distance of 2*eps of the vertex B. Note that the median from A
/// (the segment connecting A to the midpoint of BC) passes through S, since
/// this is the shortest path connecting the two endpoints. On the other
/// hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
/// the surface is a much more reasonable interpretation of the "center" of
/// this triangle.
pub fn planar_centroid(a: &Point, b: &Point, c: &Point) -> Point {
    Point((&(&a.0 + &b.0) + &c.0) * (1. / 3.))
}

impl Point {
    /// chordangle constructs a ChordAngle corresponding to the distance
    /// between the two given points. The points must be unit length.
    pub fn chordangle(&self, other: &Point) -> ChordAngle {
        ChordAngle(4f64.min((&self.0 - &other.0).norm2()))
    }
}

/*
// regularPoints generates a slice of points shaped as a regular polygon with
// the numVertices vertices, all located on a circle of the specified angular radius
// around the center. The radius is the actual distance from center to each vertex.
func regularPoints(center Point, radius s1.Angle, numVertices int) []Point {
	return regularPointsForFrame(getFrame(center), radius, numVertices)
}

// regularPointsForFrame generates a slice of points shaped as a regular polygon
// with numVertices vertices, all on a circle of the specified angular radius around
// the center. The radius is the actual distance from the center to each vertex.
func regularPointsForFrame(frame matrix3x3, radius s1.Angle, numVertices int) []Point {
	// We construct the loop in the given frame coordinates, with the center at
	// (0, 0, 1). For a loop of radius r, the loop vertices have the form
	// (x, y, z) where x^2 + y^2 = sin(r) and z = cos(r). The distance on the
	// sphere (arc length) from each vertex to the center is acos(cos(r)) = r.
	z := math.Cos(radius.Radians())
	r := math.Sin(radius.Radians())
	radianStep := 2 * math.Pi / float64(numVertices)
	var vertices []Point

	for i := 0; i < numVertices; i++ {
		angle := float64(i) * radianStep
		p := Point{r3.Vector{r * math.Cos(angle), r * math.Sin(angle), z}}
		vertices = append(vertices, Point{fromFrame(frame, p).Normalize()})
	}

	return vertices
}
*/

impl Region for Point {
    /// cap_bound returns a bounding cap for this point.
    fn cap_bound(&self) -> Cap {
        Cap::from(self)
    }

    /// rect_bound returns a bounding latitude-longitude rectangle from this point.
    fn rect_bound(&self) -> Rect {
        Rect::from(LatLng::from(self))
    }

    /// contains_cell returns false as Points do not contain any other S2 types.
    fn contains_cell(&self, _: &Cell) -> bool {
        false
    }

    /// intersects_cell reports whether this Point intersects the given cell.
    fn intersects_cell(&self, c: &Cell) -> bool {
        c.contains_point(self)
    }
}

impl Point {
    pub fn contains(&self, other: &Point) -> bool {
        self == other
    }
}

// TODO: Differences from C++
// Rotate
// Angle
// TurnAngle
// SignedArea

#[cfg(test)]
mod tests {
    use super::*;

    use std::f64::consts::PI;
    use consts::*;
    use s2::stuv::st_to_uv;

    #[test]
    fn test_origin_point() {
        assert!((Point::origin().norm() - 1.).abs() <= EPSILON);

        // The point chosen below is about 66km from the north pole towards the East
        // Siberian Sea. The purpose of the stToUV(2/3) calculation is to keep the
        // origin as far away as possible from the longitudinal edges of large
        // Cells. (The line of longitude through the chosen point is always 1/3
        // or 2/3 of the way across any Cell with longitudinal edges that it
        // passes through.)
        let p = Point::from_coords(-0.01, 0.01 * st_to_uv(2. / 3.), 1.);
        assert!(p.approx_eq(&Point::origin()));

        // Check that the origin is not too close to either pole.
        // The Earth's mean radius in kilometers (according to NASA).
        const EARTH_RADIUS_KM: f64 = 6371.01;
        assert!(Point::origin().0.z.acos() * EARTH_RADIUS_KM > 50.);
    }

    fn test_point_cross_case(expected: f64, v1: Vector, v2: Vector) {
        let p1 = Point(v1);
        let p2 = Point(v2);

        let result = p1.cross(&p2);
        assert!(f64_eq(expected, result.norm()),
                "{} != {}",
                expected,
                result.norm());
        assert!(f64_eq(0., result.0.dot(&p1.0)));
        assert!(f64_eq(0., result.0.dot(&p2.0)));
    }

    #[test]
    fn test_point_cross() {
        test_point_cross_case(1., Vector::xyz(1., 0., 0.), Vector::xyz(1., 0., 0.));
        test_point_cross_case(2., Vector::xyz(1., 0., 0.), Vector::xyz(0., 1., 0.));
        test_point_cross_case(2., Vector::xyz(0., 1., 0.), Vector::xyz(1., 0., 0.));
        test_point_cross_case(2. * 934f64.sqrt(),
                              Vector::xyz(1., 2., 3.),
                              Vector::xyz(-4., 5., -6.));
    }

    fn test_point_distance_case(expected: f64, v1: Vector, v2: Vector) {
        let p1 = Point(v1);
        let p2 = Point(v2);

        assert!(f64_eq(expected, p1.distance(&p2).0));
        assert!(f64_eq(expected, p2.distance(&p1).0));
    }

    #[test]
    fn test_point_distance() {
        test_point_distance_case(0., Vector::xyz(1., 0., 0.), Vector::xyz(1., 0., 0.));
        test_point_distance_case(PI / 2., Vector::xyz(1., 0., 0.), Vector::xyz(0., 1., 0.));
        test_point_distance_case(PI / 2., Vector::xyz(1., 0., 0.), Vector::xyz(0., 1., 1.));
        test_point_distance_case(1.2055891055045298,
                                 Vector::xyz(1., 2., 3.),
                                 Vector::xyz(2., 3., -1.));
    }
}

/*

func TestChordAngleBetweenPoints(t *testing.T) {
	for iter := 0; iter < 10; iter++ {
		m := randomFrame()
		x := m.col(0)
		y := m.col(1)
		z := m.col(2)

		if got := ChordAngleBetweenPoints(z, z).Angle(); got != 0 {
			t.Errorf("ChordAngleBetweenPoints(%v, %v) = %v, want 0", z, z, got)
		}
        if got, want := ChordAngleBetweenPoints(Point{z.Mul(-1)}, z).Angle().Radians(), math.Pi;
        !float64Near(got, want, 1e-7) {
			t.Errorf("ChordAngleBetweenPoints(%v, %v) = %v, want %v", z.Mul(-1), z, got, want)
		}
		if got, want := ChordAngleBetweenPoints(x, z).Angle().Radians(), math.Pi/2; !float64Eq(got, want) {
			t.Errorf("ChordAngleBetweenPoints(%v, %v) = %v, want %v", x, z, got, want)
		}
		w := Point{y.Add(z.Vector).Normalize()}
		if got, want := ChordAngleBetweenPoints(w, z).Angle().Radians(), math.Pi/4; !float64Eq(got, want) {
			t.Errorf("ChordAngleBetweenPoints(%v, %v) = %v, want %v", w, z, got, want)
		}
	}
}

func TestPointApproxEqual(t *testing.T) {
	tests := []struct {
		x1, y1, z1 float64
		x2, y2, z2 float64
		want       bool
	}{
		{1, 0, 0, 1, 0, 0, true},
		{1, 0, 0, 0, 1, 0, false},
		{1, 0, 0, 0, 1, 1, false},
		{1, 0, 0, -1, 0, 0, false},
		{1, 2, 3, 2, 3, -1, false},
		{1, 0, 0, 1 * (1 + epsilon), 0, 0, true},
		{1, 0, 0, 1 * (1 - epsilon), 0, 0, true},
		{1, 0, 0, 1 + epsilon, 0, 0, true},
		{1, 0, 0, 1 - epsilon, 0, 0, true},
		{1, 0, 0, 1, epsilon, 0, true},
		{1, 0, 0, 1, epsilon, epsilon, false},
		{1, epsilon, 0, 1, -epsilon, epsilon, false},
	}
	for _, test := range tests {
		p1 := Point{r3.Vector{test.x1, test.y1, test.z1}}
		p2 := Point{r3.Vector{test.x2, test.y2, test.z2}}
		if got := p1.ApproxEqual(p2); got != test.want {
			t.Errorf("%v.ApproxEqual(%v), got %v want %v", p1, p2, got, test.want)
		}
	}
}

var (
	pz   = Point{r3.Vector{0, 0, 1}}
	p000 = Point{r3.Vector{1, 0, 0}}
	p045 = Point{r3.Vector{1, 1, 0}}
	p090 = Point{r3.Vector{0, 1, 0}}
	p180 = Point{r3.Vector{-1, 0, 0}}
	// Degenerate triangles.
	pr = Point{r3.Vector{0.257, -0.5723, 0.112}}
	pq = Point{r3.Vector{-0.747, 0.401, 0.2235}}

	// For testing the Girard area fall through case.
	g1 = Point{r3.Vector{1, 1, 1}}
	g2 = Point{g1.Add(pr.Mul(1e-15)).Normalize()}
	g3 = Point{g1.Add(pq.Mul(1e-15)).Normalize()}
)

func TestPointArea(t *testing.T) {
	epsilon := 1e-10
	tests := []struct {
		a, b, c  Point
		want     float64
		nearness float64
	}{
		{p000, p090, pz, math.Pi / 2.0, 0},
		// This test case should give 0 as the epsilon, but either Go or C++'s value for Pi,
		// or the accuracy of the multiplications along the way, cause a difference ~15 decimal
		// places into the result, so it is not quite a difference of 0.
		{p045, pz, p180, 3.0 * math.Pi / 4.0, 1e-14},
		// Make sure that Area has good *relative* accuracy even for very small areas.
		{Point{r3.Vector{epsilon, 0, 1}}, Point{r3.Vector{0, epsilon, 1}}, pz, 0.5 * epsilon * epsilon, 1e-14},
		// Make sure that it can handle degenerate triangles.
		{pr, pr, pr, 0.0, 0},
		{pr, pq, pr, 0.0, 1e-15},
		{p000, p045, p090, 0.0, 0},
		// Try a very long and skinny triangle.
		{p000, Point{r3.Vector{1, 1, epsilon}}, p090, 5.8578643762690495119753e-11, 1e-9},
		// TODO(roberts):
		// C++ includes a 10,000 loop of perterbations to test out the Girard area
		// computation is less than some noise threshold.
		// Do we need that many? Will one or two suffice?
		{g1, g2, g3, 0.0, 1e-15},
	}
	for _, test := range tests {
		if got := PointArea(test.a, test.b, test.c); !float64Near(got, test.want, test.nearness) {
			t.Errorf("PointArea(%v, %v, %v), got %v want %v", test.a, test.b, test.c, got, test.want)
		}
	}
}

func TestPointAreaQuarterHemisphere(t *testing.T) {
	tests := []struct {
		a, b, c, d, e Point
		want          float64
	}{
		// Triangles with near-180 degree edges that sum to a quarter-sphere.
		{Point{r3.Vector{1, 0.1 * epsilon, epsilon}}, p000, p045, p180, pz, math.Pi},
		// Four other triangles that sum to a quarter-sphere.
		{Point{r3.Vector{1, 1, epsilon}}, p000, p045, p180, pz, math.Pi},
		// TODO(roberts):
		// C++ Includes a loop of 100 perturbations on a hemisphere for more tests.
	}
	for _, test := range tests {
		area := PointArea(test.a, test.b, test.c) +
			PointArea(test.a, test.c, test.d) +
			PointArea(test.a, test.d, test.e) +
			PointArea(test.a, test.e, test.b)

		if !float64Eq(area, test.want) {
			t.Errorf("Adding up 4 quarter hemispheres with PointArea(), got %v want %v", area, test.want)
		}
	}
}

func TestPointPlanarCentroid(t *testing.T) {
	tests := []struct {
		name             string
		p0, p1, p2, want Point
	}{
		{
			name: "xyz axis",
			p0:   Point{r3.Vector{0, 0, 1}},
			p1:   Point{r3.Vector{0, 1, 0}},
			p2:   Point{r3.Vector{1, 0, 0}},
			want: Point{r3.Vector{1. / 3, 1. / 3, 1. / 3}},
		},
		{
			name: "Same point",
			p0:   Point{r3.Vector{1, 0, 0}},
			p1:   Point{r3.Vector{1, 0, 0}},
			p2:   Point{r3.Vector{1, 0, 0}},
			want: Point{r3.Vector{1, 0, 0}},
		},
	}

	for _, test := range tests {
		got := PlanarCentroid(test.p0, test.p1, test.p2)
		if !got.ApproxEqual(test.want) {
			t.Errorf("%s: PlanarCentroid(%v, %v, %v) = %v, want %v", test.name, test.p0, test.p1, test.p2, got, test.want)
		}
	}
}

func TestPointTrueCentroid(t *testing.T) {
	// Test TrueCentroid with very small triangles. This test assumes that
	// the triangle is small enough so that it is nearly planar.
	// The centroid of a planar triangle is at the intersection of its
	// medians, which is two-thirds of the way along each median.
	for i := 0; i < 100; i++ {
		f := randomFrame()
		p := f.col(0)
		x := f.col(1)
		y := f.col(2)
		d := 1e-4 * math.Pow(1e-4, randomFloat64())

		// Make a triangle with two equal sides.
		p0 := Point{p.Sub(x.Mul(d)).Normalize()}
		p1 := Point{p.Add(x.Mul(d)).Normalize()}
		p2 := Point{p.Add(y.Mul(d * 3)).Normalize()}
		want := Point{p.Add(y.Mul(d)).Normalize()}

		got := TrueCentroid(p0, p1, p2).Normalize()
		if got.Distance(want.Vector) >= 2e-8 {
			t.Errorf("TrueCentroid(%v, %v, %v).Normalize() = %v, want %v", p0, p1, p2, got, want)
		}

		// Make a triangle with a right angle.
		p0 = p
		p1 = Point{p.Add(x.Mul(d * 3)).Normalize()}
		p2 = Point{p.Add(y.Mul(d * 6)).Normalize()}
		want = Point{p.Add(x.Add(y.Mul(2)).Mul(d)).Normalize()}

		got = TrueCentroid(p0, p1, p2).Normalize()
		if got.Distance(want.Vector) >= 2e-8 {
			t.Errorf("TrueCentroid(%v, %v, %v).Normalize() = %v, want %v", p0, p1, p2, got, want)
		}
	}
}

func TestPointRegularPoints(t *testing.T) {
	// Conversion to/from degrees has a little more variability than the default epsilon.
	const epsilon = 1e-13
	center := PointFromLatLng(LatLngFromDegrees(80, 135))
	radius := s1.Degree * 20
	pts := regularPoints(center, radius, 4)

	if len(pts) != 4 {
		t.Errorf("regularPoints with 4 vertices should have 4 vertices, got %d", len(pts))
	}

	lls := []LatLng{
		LatLngFromPoint(pts[0]),
		LatLngFromPoint(pts[1]),
		LatLngFromPoint(pts[2]),
		LatLngFromPoint(pts[3]),
	}
	cll := LatLngFromPoint(center)

	// Make sure that the radius is correct.
	wantDist := 20.0
	for i, ll := range lls {
		if got := cll.Distance(ll).Degrees(); !float64Near(got, wantDist, epsilon) {
			t.Errorf("Vertex %d distance from center = %v, want %v", i, got, wantDist)
		}
	}

	// Make sure the angle between each point is correct.
	wantAngle := math.Pi / 2
	for i := 0; i < len(pts); i++ {
		// Mod the index by 4 to wrap the values at each end.
		v0, v1, v2 := pts[(4+i+1)%4], pts[(4+i)%4], pts[(4+i-1)%4]
		if got := float64(v0.Sub(v1.Vector).Angle(v2.Sub(v1.Vector))); !float64Eq(got, wantAngle) {
			t.Errorf("(%v-%v).Angle(%v-%v) = %v, want %v", v0, v1, v1, v2, got, wantAngle)
		}
	}

	// Make sure that all edges of the polygon have the same length.
	wantLength := 27.990890717782829
	for i := 0; i < len(lls); i++ {
		ll1, ll2 := lls[i], lls[(i+1)%4]
		if got := ll1.Distance(ll2).Degrees(); !float64Near(got, wantLength, epsilon) {
			t.Errorf("%v.Distance(%v) = %v, want %v", ll1, ll2, got, wantLength)
		}
	}

	// Spot check an actual coordinate now that we know the points are spaced
	// evenly apart at the same angles and radii.
	if got, want := lls[0].Lat.Degrees(), 62.162880741097204; !float64Near(got, want, epsilon) {
		t.Errorf("%v.Lat = %v, want %v", lls[0], got, want)
	}
	if got, want := lls[0].Lng.Degrees(), 103.11051028343407; !float64Near(got, want, epsilon) {
		t.Errorf("%v.Lng = %v, want %v", lls[0], got, want)
	}
}

func TestPointRegion(t *testing.T) {
	p := Point{r3.Vector{1, 0, 0}}
	r := Point{r3.Vector{1, 0, 0}}
	if !r.Contains(p) {
		t.Errorf("%v.Contains(%v) = false, want true", r, p)
	}
	if !r.Contains(r) {
		t.Errorf("%v.Contains(%v) = false, want true", r, r)
	}
	if s := (Point{r3.Vector{1, 0, 1}}); r.Contains(s) {
		t.Errorf("%v.Contains(%v) = true, want false", r, s)
	}
	if got, want := r.CapBound(), CapFromPoint(p); !got.ApproxEqual(want) {
		t.Errorf("%v.CapBound() = %v, want %v", r, got, want)
	}
	if got, want := r.RectBound(), RectFromLatLng(LatLngFromPoint(p)); !rectsApproxEqual(got, want, epsilon, epsilon) {
		t.Errorf("%v.RectBound() = %v, want %v", r, got, want)
	}

	// The leaf cell containing a point is still much larger than the point.
	cell := CellFromPoint(p)
	if r.ContainsCell(cell) {
		t.Errorf("%v.ContainsCell(%v) = true, want false", r, cell)
	}
	if !r.IntersectsCell(cell) {
		t.Errorf("%v.IntersectsCell(%v) = false, want true", r, cell)
	}
}

func BenchmarkPointArea(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PointArea(p000, p090, pz)
	}
}

func BenchmarkPointAreaGirardCase(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PointArea(g1, g2, g3)
	}
}
*/
