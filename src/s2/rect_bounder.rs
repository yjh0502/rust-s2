// Copyright 2017 Google Inc. All rights reserved.
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

// RectBounder is used to compute a bounding rectangle that contains all edges
// defined by a vertex chain (v0, v1, v2, ...). All vertices must be unit length.
// Note that the bounding rectangle of an edge can be larger than the bounding
// rectangle of its endpoints, e.g. consider an edge that passes through the North Pole.
//
// The bounds are calculated conservatively to account for numerical errors
// when points are converted to LatLngs. More precisely, this function
// guarantees the following:
// Let L be a closed edge chain (Loop) such that the interior of the loop does
// not contain either pole. Now if P is any point such that L.ContainsPoint(P),
// then RectBounder(L).ContainsPoint(LatLngFromPoint(P)).

use std::f64::consts::{FRAC_PI_2, PI};

use crate::consts::DBL_EPSILON;
use crate::r1;
use crate::r3::vector::Vector;
use crate::s1;
use crate::s2::latlng::LatLng;
use crate::s2::point::Point;
use crate::s2::rect::Rect;

static Z_AXIS_POINT: Point = Point(Vector {
    x: 0.,
    y: 0.,
    z: 1.,
});

// A vertex in a Point chain.
struct Vertex {
    point: Point,
    ll: LatLng,
}

// RectBounder is used to compute a bounding rectangle that contains all edges
// defined by a vertex chain (v0, v1, v2, ...). All vertices must be unit length.
// Note that the bounding rectangle of an edge can be larger than the bounding
// rectangle of its endpoints, e.g. consider an edge that passes through the North Pole.
//
// The bounds are calculated conservatively to account for numerical errors
// when points are converted to LatLngs. More precisely, this function
// guarantees the following:
// Let L be a closed edge chain (Loop) such that the interior of the loop does
// not contain either pole. Now if P is any point such that L.ContainsPoint(P),
// then RectBounder(L).ContainsPoint(LatLngFromPoint(P)).
pub struct RectBounder {
    // The previous vertex in the chain.
    a: Option<Vertex>,

    // The current bounding rectangle.
    bound: Rect,
}

impl Default for RectBounder {
    fn default() -> Self {
        Self::new()
    }
}

impl RectBounder {
    pub fn new() -> Self {
        RectBounder {
            a: None,
            bound: Rect::empty(),
        }
    }

    // Adds the given point to the chain. The Point must be unit length.
    pub fn add_point(&mut self, b: &Point) {
        let b_ll = LatLng::from(b);

        match &self.a {
            None => {
                self.bound = &self.bound + &b_ll;
            }
            Some(a) => {
                // First compute the cross product N = A x B robustly. This is the normal
                // to the great circle through A and B. We don't use RobustSign
                // since that method returns an arbitrary vector orthogonal to A if the two
                // vectors are proportional, and we want the zero vector in that case.
                let n = Point((a.point - *b).0.cross(&(a.point + *b).0)); // N = 2 * (A x B)

                // The relative error in N gets large as its norm gets very small (i.e.,
                // when the two points are nearly identical or antipodal). We handle this
                // by choosing a maximum allowable error, and if the error is greater than
                // this we fall back to a different technique. Since it turns out that
                // the other sources of error in converting the normal to a maximum
                // latitude add up to at most 1.16 * dblEpsilon, and it is desirable to
                // have the total error be a multiple of dblEpsilon, we have chosen to
                // limit the maximum error in the normal to be 3.84 * dblEpsilon.
                // It is possible to show that the error is less than this when
                //
                // n.Norm() >= 8 * sqrt(3) / (3.84 - 0.5 - sqrt(3)) * dblEpsilon
                //          = 1.91346e-15 (about 8.618 * dblEpsilon)
                let n_norm = n.norm();
                if n_norm < 1.91346e-15 {
                    // A and B are either nearly identical or nearly antipodal (to within
                    // 4.309 * dblEpsilon, or about 6 nanometers on the earth's surface).
                    if a.point.0.dot(&b.0) < 0. {
                        // The two points are nearly antipodal. The easiest solution is to
                        // assume that the edge between A and B could go in any direction
                        // around the sphere.
                        self.bound = Rect::full();
                    } else {
                        // The two points are nearly identical (to within 4.309 * dblEpsilon).
                        // In this case we can just use the bounding rectangle of the points,
                        // since after the expansion done by GetBound this Rect is
                        // guaranteed to include the (lat,lng) values of all points along AB.
                        self.bound = self.bound.union(&Rect::from_point_pair(&a.ll, &b_ll));
                    }
                } else {
                    // Compute the longitude range spanned by AB.
                    let mut lng_ab =
                        s1::interval::Interval::from_point_pair(a.ll.lng.rad(), b_ll.lng.rad());
                    if lng_ab.len() >= PI - 2. * DBL_EPSILON {
                        // The points lie on nearly opposite lines of longitude to within the
                        // maximum error of the calculation. The easiest solution is to assume
                        // that AB could go on either side of the pole.
                        lng_ab = s1::interval::FULL;
                    }

                    // Next we compute the latitude range spanned by the edge AB.  We start
                    // with the range spanning the two endpoints of the edge:
                    let mut lat_ab =
                        r1::interval::Interval::from_point_pair(a.ll.lat.rad(), b_ll.lat.rad());

                    // This is the desired range unless the edge AB crosses the plane
                    // through N and the Z-axis (which is where the great circle through A
                    // and B attains its minimum and maximum latitudes).  To test whether AB
                    // crosses this plane, we compute a vector M perpendicular to this
                    // plane and then project A and B onto it.
                    let m = Point(n.0.cross(&Z_AXIS_POINT.0));
                    let m_a = m.0.dot(&a.point.0);
                    let m_b = m.0.dot(&b.0);

                    // We want to test the signs of "mA" and "mB", so we need to bound
                    // the error in these calculations. It is possible to show that the
                    // total error is bounded by
                    //
                    // (1 + sqrt(3)) * dblEpsilon * nNorm + 8 * sqrt(3) * (dblEpsilon**2)
                    //   = 6.06638e-16 * nNorm + 6.83174e-31

                    let m_error = 6.06638e-16 * n_norm + 6.83174e-31;
                    if m_a * m_b < 0. || m_a.abs() <= m_error || m_b.abs() <= m_error {
                        // Minimum/maximum latitude *may* occur in the edge interior.
                        //
                        // The maximum latitude is 90 degrees minus the latitude of N. We
                        // compute this directly using atan2 in order to get maximum accuracy
                        // near the poles.
                        //
                        // Our goal is compute a bound that contains the computed latitudes of
                        // all Points P that pass the point-in-polygon containment test.
                        // There are three sources of error we need to consider:
                        // - the directional error in N (at most 3.84 * dblEpsilon)
                        // - converting N to a maximum latitude
                        // - computing the latitude of the test point P
                        // The latter two sources of error are at most 0.955 * dblEpsilon
                        // individually, but it is possible to show by a more complex analysis
                        // that together they can add up to at most 1.16 * dblEpsilon, for a
                        // total error of 5 * dblEpsilon.
                        //
                        // We add 3 * dblEpsilon to the bound here, and GetBound() will pad
                        // the bound by another 2 * dblEpsilon.
                        let max_lat = ((n.0.x * n.0.x + n.0.y * n.0.y).sqrt().atan2(n.0.z.abs())
                            + 3. * DBL_EPSILON)
                            .min(FRAC_PI_2);

                        // In order to get tight bounds when the two points are close together,
                        // we also bound the min/max latitude relative to the latitudes of the
                        // endpoints A and B. First we compute the distance between A and B,
                        // and then we compute the maximum change in latitude between any two
                        // points along the great circle that are separated by this distance.
                        // This gives us a latitude change "budget". Some of this budget must
                        // be spent getting from A to B; the remainder bounds the round-trip
                        // distance (in latitude) from A or B to the min or max latitude
                        // attained along the edge AB.
                        let lat_budget_z = 0.5 * (a.point.0 - b.0).norm() * max_lat.sin();
                        let lat_budget =
                            2. * (((1. + 4. * DBL_EPSILON) * lat_budget_z).min(1.0)).asin();

                        let max_delta = 0.5 * (lat_budget - lat_ab.len()) + DBL_EPSILON;

                        // Test whether AB passes through the point of maximum latitude or
                        // minimum latitude. If the dot product(s) are small enough then the
                        // result may be ambiguous.
                        if m_a <= m_error && m_b >= -m_error {
                            lat_ab.hi = max_lat.min(lat_ab.hi + max_delta);
                        }
                        if m_b <= m_error && m_a >= -m_error {
                            lat_ab.lo = (-max_lat).max(lat_ab.lo - max_delta);
                        }
                    }
                    self.bound = self.bound.union(&Rect {
                        lat: lat_ab,
                        lng: lng_ab,
                    });
                }
            }
        }

        // Update the previous vertex.
        self.a = Some(Vertex {
            point: *b,
            ll: b_ll,
        });
    }

    // Returns the bounding rectangle of the edge chain that connects the
    // vertices defined so far.  This bound satisfies the guarantee made
    // above, i.e. if the edge chain defines a loop, then the bound contains
    // the LatLng coordinates of all Points contained by the loop.
    pub fn get_bound(&self) -> Rect {
        // To save time, we ignore numerical errors in the computed LatLngs while
        // accumulating the bounds and then account for them here.
        //
        // LatLng(Point) has a maximum error of 0.955 * DBL_EPSILON in latitude.
        // In the worst case, we might have rounded "inwards" when computing the
        // bound and "outwards" when computing the latitude of a contained point P,
        // therefore we expand the latitude bounds by 2 * DBL_EPSILON in each
        // direction.  (A more complex analysis shows that 1.5 * DBL_EPSILON is
        // enough, but the expansion amount should be a multiple of DBL_EPSILON in
        // order to avoid rounding errors during the expansion itself.)
        //
        // LatLng(Point) has a maximum error of DBL_EPSILON in longitude, which
        // is simply the maximum rounding error for results in the range [-Pi, Pi].
        // This is true because the Gnu implementation of atan2() comes from the IBM
        // Accurate Mathematical Library, which implements correct rounding for this
        // intrinsic (i.e., it returns the infinite precision result rounded to the
        // nearest representable value, with ties rounded to even values).  This
        // implies that we don't need to expand the longitude bounds at all, since
        // we only guarantee that the bound contains the *rounded* latitudes of
        // contained points.  The *true* latitudes of contained points may lie up to
        // DBL_EPSILON outside of the returned bound.
        let expansion = LatLng::from_radians(2. * DBL_EPSILON, 0.);
        self.bound.expanded(&expansion).polar_closure()
    }
}

// Expands a bound returned by GetBound() so that it is guaranteed to
// contain the bounds of any subregion whose bounds are computed using
// this class.  For example, consider a loop L that defines a square.
// get_bound() ensures that if a point P is contained by this square, then
// LatLng(P) is contained by the bound.  But now consider a diamond
// shaped loop S contained by L.  It is possible that get_bound() returns a
// *larger* bound for S than it does for L, due to rounding errors.  This
// method expands the bound for L so that it is guaranteed to contain the
// bounds of any subregion S.
//
// More precisely, if L is a loop that does not contain either pole, and S
// is a loop such that L.Contains(S), then
//
//   ExpandForSubregions(RectBounder(L)).Contains(RectBounder(S)).
pub fn expand_for_subregions(bound: &Rect) -> Rect {
    // Empty bounds don't need expansion.
    if bound.is_empty() {
        return bound.clone();
    }

    // First we need to check whether the bound B contains any nearly-antipodal
    // points (to within 4.309 * dblEpsilon). If so then we need to return
    // FullRect, since the subregion might have an edge between two
    // such points, and AddPoint returns Full for such edges. Note that
    // this can happen even if B is not Full for example, consider a loop
    // that defines a 10km strip straddling the equator extending from
    // longitudes -100 to +100 degrees.
    //
    // It is easy to check whether B contains any antipodal points, but checking
    // for nearly-antipodal points is trickier. Essentially we consider the
    // original bound B and its reflection through the origin B', and then test
    // whether the minimum distance between B and B' is less than 4.309 * dblEpsilon.

    // lngGap is a lower bound on the longitudinal distance between B and its
    // reflection B'. (2.5 * dblEpsilon is the maximum combined error of the
    // endpoint longitude calculations and the Length call.)
    let lng_gap = 0f64.max(PI - bound.lng.len() - 2.5 * DBL_EPSILON);

    // minAbsLat is the minimum distance from B to the equator (if zero or
    // negative, then B straddles the equator).
    let min_abs_lat = bound.lat.lo.max(-bound.lat.hi);

    // latGapSouth and latGapNorth measure the minimum distance from B to the
    // south and north poles respectively.
    let lat_gap_south = FRAC_PI_2 + bound.lat.lo;
    let lat_gap_north = FRAC_PI_2 - bound.lat.hi;

    if min_abs_lat >= 0. {
        // The bound B does not straddle the equator. In this case the minimum
        // distance is between one endpoint of the latitude edge in B closest to
        // the equator and the other endpoint of that edge in B'. The latitude
        // distance between these two points is 2*minAbsLat, and the longitude
        // distance is lngGap. We could compute the distance exactly using the
        // Haversine formula, but then we would need to bound the errors in that
        // calculation. Since we only need accuracy when the distance is very
        // small (close to 4.309 * dblEpsilon), we substitute the Euclidean
        // distance instead. This gives us a right triangle XYZ with two edges of
        // length x = 2*minAbsLat and y ~= lngGap. The desired distance is the
        // length of the third edge z, and we have
        //
        //         z  ~=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
        //
        // Therefore the region may contain nearly antipodal points only if
        //
        //  2*minAbsLat + lngGap  <  sqrt(2) * 4.309 * dblEpsilon
        //                        ~= 1.354e-15
        //
        // Note that because the given bound B is conservative, minAbsLat and
        // lngGap are both lower bounds on their true values so we do not need
        // to make any adjustments for their errors.
        if 2. * min_abs_lat + lng_gap < 1.354e-15 {
            return Rect::full();
        }
    } else if lng_gap >= FRAC_PI_2 {
        // B spans at most Pi/2 in longitude. The minimum distance is always
        // between one corner of B and the diagonally opposite corner of B'. We
        // use the same distance approximation that we used above; in this case
        // we have an obtuse triangle XYZ with two edges of length x = latGapSouth
        // and y = latGapNorth, and angle Z >= Pi/2 between them. We then have
        //
        //         z  >=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
        //
        // Unlike the case above, latGapSouth and latGapNorth are not lower bounds
        // (because of the extra addition operation, and because math.Pi/2 is not
        // exactly equal to Pi/2); they can exceed their true values by up to
        // 0.75 * dblEpsilon. Putting this all together, the region may contain
        // nearly antipodal points only if
        //
        //   latGapSouth + latGapNorth  <  (sqrt(2) * 4.309 + 1.5) * dblEpsilon
        //                              ~= 1.687e-15
        if lat_gap_south + lat_gap_north < 1.687e-15 {
            return Rect::full();
        }
    } else {
        // Otherwise we know that (1) the bound straddles the equator and (2) its
        // width in longitude is at least Pi/2. In this case the minimum
        // distance can occur either between a corner of B and the diagonally
        // opposite corner of B' (as in the case above), or between a corner of B
        // and the opposite longitudinal edge reflected in B'. It is sufficient
        // to only consider the corner-edge case, since this distance is also a
        // lower bound on the corner-corner distance when that case applies.

        // Consider the spherical triangle XYZ where X is a corner of B with
        // minimum absolute latitude, Y is the closest pole to X, and Z is the
        // point closest to X on the opposite longitudinal edge of B'. This is a
        // right triangle (Z = Pi/2), and from the spherical law of sines we have
        //
        //     sin(z) / sin(Z)  =  sin(y) / sin(Y)
        //     sin(maxLatGap) / 1  =  sin(dMin) / sin(lngGap)
        //     sin(dMin)  =  sin(maxLatGap) * sin(lngGap)
        //
        // where "maxLatGap" = max(latGapSouth, latGapNorth) and "dMin" is the
        // desired minimum distance. Now using the facts that sin(t) >= (2/Pi)*t
        // for 0 <= t <= Pi/2, that we only need an accurate approximation when
        // at least one of "maxLatGap" or lngGap is extremely small (in which
        // case sin(t) ~= t), and recalling that "maxLatGap" has an error of up
        // to 0.75 * dblEpsilon, we want to test whether
        //
        //   maxLatGap * lngGap  <  (4.309 + 0.75) * (Pi/2) * dblEpsilon
        //                       ~= 1.765e-15
        if lat_gap_south.max(lat_gap_north) * lng_gap < 1.765e-15 {
            return Rect::full();
        }
    }

    // Next we need to check whether the subregion might contain any edges that
    // span (math.Pi - 2 * dblEpsilon) radians or more in longitude, since AddPoint
    // sets the longitude bound to Full in that case. This corresponds to
    // testing whether (lngGap <= 0) in lngExpansion below.

    // Otherwise, the maximum latitude error in AddPoint is 4.8 * dblEpsilon.
    // In the worst case, the errors when computing the latitude bound for a
    // subregion could go in the opposite direction as the errors when computing
    // the bound for the original region, so we need to double this value.
    // (More analysis shows that it's okay to round down to a multiple of
    // dblEpsilon.)
    //
    // For longitude, we rely on the fact that atan2 is correctly rounded and
    // therefore no additional bounds expansion is necessary.
    let lat_expansion = 9. * DBL_EPSILON;
    let lng_expansion = if lng_gap <= 0. { PI } else { 0. };
    bound
        .expanded(&LatLng::from_radians(lat_expansion, lng_expansion))
        .polar_closure()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::consts::f64_near;
    use crate::r3::vector::Vector;
    use crate::s2::random;
    use crate::s2::rect::VALID_RECT_LAT_RANGE;
    use rand::Rng;

    pub fn max_error_for_tests() -> LatLng {
        // The maximum error in the latitude calculation is
        //    3.84 * DBL_EPSILON   for the cross product calculation (see above)
        //    0.96 * DBL_EPSILON   for the Latitude() calculation
        //    5    * DBL_EPSILON   added by add_point/get_bound to compensate for error
        //    ------------------
        //    9.80 * DBL_EPSILON   maximum error in result
        //
        // The maximum error in the longitude calculation is DBL_EPSILON. get_bound()
        // does not do any expansion because this isn't necessary in order to
        // bound the *rounded* longitudes of contained points.
        return LatLng::from_radians(10. * DBL_EPSILON, 1. * DBL_EPSILON);
    }

    fn rect_bound_for_points(a: &Point, b: &Point) -> Rect {
        let mut bounder = RectBounder::new();
        bounder.add_point(a);
        bounder.add_point(b);
        bounder.get_bound()
    }

    #[test]
    fn test_max_latitude_simple() {
        let cube_lat = (1. / 3.0f64.sqrt()).asin(); // 35.26 degrees
        let cube_lat_rect = Rect {
            lat: r1::interval::Interval::new(-cube_lat, cube_lat),
            lng: s1::interval::Interval::new(-PI / 4., PI / 4.),
        };

        fn test(a: &Point, b: &Point, expected: &Rect) {
            let got = rect_bound_for_points(a, b);
            assert!(
                got.approx_eq(expected),
                "RectBounder for points ({}, {}) near max lat failed: got {}, want {}",
                a,
                b,
                got,
                expected
            );
        }

        // Check cases where the min/max latitude is attained at a vertex.
        test(
            &Point(Vector {
                x: 1.,
                y: 1.,
                z: 1.,
            }),
            &Point(Vector {
                x: 1.,
                y: -1.,
                z: -1.,
            }),
            &cube_lat_rect,
        );
        test(
            &Point(Vector {
                x: 1.,
                y: -1.,
                z: 1.,
            }),
            &Point(Vector {
                x: 1.,
                y: 1.,
                z: -1.,
            }),
            &cube_lat_rect,
        );
    }

    #[test]
    fn test_rect_bounder_max_latitude_edge_interior() {
        // Check cases where the min/max latitude occurs in the edge interior.
        // These tests expect the result to be pretty close to the middle of the
        // allowable error range (i.e., by adding 0.5 * kRectError).

        // Max latitude, CW edge
        assert_eq!(
            rect_bound_for_points(
                &Point(Vector {
                    x: 1.,
                    y: 1.,
                    z: 1.
                }),
                &Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.
                })
            )
            .lat
            .hi,
            PI / 4. + 0.5 * max_error_for_tests().lat.rad()
        );

        // Min latitude, CW edge
        assert_eq!(
            rect_bound_for_points(
                &Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: -1.
                }),
                &Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: -1.
                })
            )
            .lat
            .lo,
            -PI / 4. - 0.5 * max_error_for_tests().lat.rad()
        );

        // Max latitude, CCW edge
        assert_eq!(
            rect_bound_for_points(
                &Point(Vector {
                    x: 1.,
                    y: -1.,
                    z: 1.
                }),
                &Point(Vector {
                    x: 1.,
                    y: 1.,
                    z: 1.
                })
            )
            .lat
            .hi,
            PI / 4. + 0.5 * max_error_for_tests().lat.rad()
        );

        // Min latitude, CCW edge
        assert_eq!(
            rect_bound_for_points(
                &Point(Vector {
                    x: -1.,
                    y: 1.,
                    z: -1.
                }),
                &Point(Vector {
                    x: -1.,
                    y: -1.,
                    z: -1.
                })
            )
            .lat
            .lo,
            -PI / 4. - 0.5 * max_error_for_tests().lat.rad()
        );

        // Check cases where the edge passes through one of the poles.
        assert_eq!(
            rect_bound_for_points(
                &Point(Vector {
                    x: 0.3,
                    y: 0.4,
                    z: 1.
                }),
                &Point(Vector {
                    x: -0.3,
                    y: -0.4,
                    z: 1.
                })
            )
            .lat
            .hi,
            FRAC_PI_2
        );
        assert_eq!(
            rect_bound_for_points(
                &Point(Vector {
                    x: 0.3,
                    y: 0.4,
                    z: -1.
                }),
                &Point(Vector {
                    x: -0.3,
                    y: -0.4,
                    z: -1.
                })
            )
            .lat
            .lo,
            -FRAC_PI_2
        );
    }

    #[test]
    fn test_rect_bounder_max_latitude_random() {
        // Check that the maximum latitude of edges is computed accurately to within
        // 3 * dblEpsilon (the expected maximum error). We concentrate on maximum
        // latitudes near the equator and north pole since these are the extremes.

        let mut rng = random::rng();
        for _ in 0..100 {
            // Construct a right-handed coordinate frame (U,V,W) such that U points
            // slightly above the equator, V points at the equator, and W is slightly
            // offset from the north pole.
            let mut u = random::point(&mut rng);
            u.0.z = DBL_EPSILON * 1e-6 * 1e12f64.powf(rng.gen::<f64>());
            u = u.normalize();

            let v = Z_AXIS_POINT.cross(&u).normalize();
            let w = u.cross(&v).normalize();

            // Construct a line segment AB that passes through U, and check that the
            // maximum latitude of this segment matches the latitude of U.
            let a = u - (v * rng.gen::<f64>()).normalize();
            let b = u + (v * rng.gen::<f64>()).normalize();
            let ab_bound = rect_bound_for_points(&a, &b);
            assert!(
                f64_near(
                    u.latitude().rad(),
                    ab_bound.lat.hi,
                    max_error_for_tests().lat.rad()
                ),
                "bound for line AB not near enough to the latitude of point {}. got {}, want {}",
                u,
                u.latitude().rad(),
                ab_bound.lat.hi
            );

            // Construct a line segment CD that passes through W, and check that the
            // maximum latitude of this segment matches the latitude of W.
            let c = w - (v * rng.gen::<f64>()).normalize();
            let d = w + (v * rng.gen::<f64>()).normalize();
            let cd_bound = rect_bound_for_points(&c, &d);
            assert!(
                f64_near(
                    w.latitude().rad(),
                    cd_bound.lat.hi,
                    max_error_for_tests().lat.rad()
                ),
                "bound for line CD not near enough to the lat of point {}. got {}, want {}",
                w,
                w.latitude().rad(),
                cd_bound.lat.hi
            );
        }
    }

    #[test]
    fn test_rect_bounder_expand_for_subregions() {
        // Test the full and empty bounds.
        assert!(
            expand_for_subregions(&Rect::full()).is_full(),
            "Subregion Bound of full rect should be full"
        );
        assert!(
            expand_for_subregions(&Rect::empty()).is_empty(),
            "Subregion Bound of empty rect should be empty"
        );

        fn get_bound(x_lat: f64, x_lng: f64, y_lat: f64, y_lng: f64) -> Rect {
            let input = Rect::from_point_pair(
                &LatLng::from_radians(x_lat, x_lng),
                &LatLng::from_radians(y_lat, y_lng),
            );
            let output = expand_for_subregions(&input);

            // Test that the bound is actually expanded.
            assert!(
                output.contains(&input),
                "Subregion bound of ({}, {}, {}, {}) should contain original rect",
                x_lat,
                x_lng,
                y_lat,
                y_lng
            );

            if input.lat == VALID_RECT_LAT_RANGE {
                assert!(
                    !input.lat.contains_interval(&output.lat),
                    "Subregion bound of ({}, {}, {}, {}) shouldn't be contained by original rect",
                    x_lat,
                    x_lng,
                    y_lat,
                    y_lng
                );
            }
            output
        }

        fn test_full(x_lat: f64, x_lng: f64, y_lat: f64, y_lng: f64, expect_full: bool) {
            let output = get_bound(x_lat, x_lng, y_lat, y_lng);

            // We check the various situations where the bound contains nearly-antipodal points. The tests are organized into pairs
            // where the two bounds are similar except that the first bound meets the nearly-antipodal criteria while the second does not.
            assert!(
                output.is_full() == expect_full,
                "Subregion Bound of ({}, {}, {}, {}).IsFull should be {}",
                x_lat,
                x_lng,
                y_lat,
                y_lng,
                expect_full
            );
        }

        // First we check the various situations where the bound contains
        // nearly-antipodal points.  The tests are organized into pairs where the
        // two bounds are similar except that the first bound meets the
        // nearly-antipodal criteria while the second does not.

        // Cases where the bound does not straddle the equator (but almost does),
        // and spans nearly 180 degrees in longitude.
        test_full(3e-16, 0., 1e-14, PI, true);
        test_full(9e-16, 0., 1e-14, PI, false);
        test_full(1e-16, 7e-16, 1e-14, PI, true);
        test_full(3e-16, 14e-16, 1e-14, PI, false);
        test_full(1e-100, 14e-16, 1e-14, PI, true);
        test_full(1e-100, 22e-16, 1e-14, PI, false);

        // Cases where the bound spans at most 90 degrees in longitude, and almost
        // 180 degrees in latitude.  Note that DBL_EPSILON is about 2.22e-16, which
        // implies that the double-precision value just below Pi/2 can be written as
        // (M_PI_2. - 2e-16).
        test_full(-FRAC_PI_2, -1e-15, FRAC_PI_2 - 7e-16, 0., true);
        test_full(-FRAC_PI_2, -1e-15, FRAC_PI_2 - 30e-16, 0., false);
        test_full(-FRAC_PI_2 + 4e-16, 0., FRAC_PI_2 - 2e-16, 1e-7, true);
        test_full(-FRAC_PI_2 + 30e-16, 0., FRAC_PI_2, 1e-7, false);
        test_full(-FRAC_PI_2 + 4e-16, 0., FRAC_PI_2 - 4e-16, FRAC_PI_2, true);
        test_full(FRAC_PI_2, 0., FRAC_PI_2 - 30e-16, FRAC_PI_2, false);

        // Cases where the bound straddles the equator and spans more than 90
        // degrees in longitude.  These are the cases where the critical distance is
        // between a corner of the bound and the opposite longitudinal edge.  Unlike
        // the cases above, here the bound may contain nearly-antipodal points (to
        // within 3.055 * DBL_EPSILON) even though the latitude and longitude ranges
        // are both significantly less than (Pi - 3.055 * DBL_EPSILON).
        test_full(-FRAC_PI_2, 0., FRAC_PI_2 - 1e-8, PI - 1e-7, true);
        test_full(-FRAC_PI_2, 0., FRAC_PI_2 - 1e-7, PI - 1e-7, false);
        test_full(-FRAC_PI_2 + 1e-12, -PI + 1e-4, FRAC_PI_2, 0., true);
        test_full(-FRAC_PI_2 + 1e-11, -PI + 1e-4, FRAC_PI_2, 0., true);

        fn test_rect(x_lat: f64, x_lng: f64, y_lat: f64, y_lng: f64, expect_rect: &Rect) {
            let output = get_bound(x_lat, x_lng, y_lat, y_lng);
            assert!(
                output.approx_eq(expect_rect),
                "Subregion Bound of ({}, {}, {}, {}) = {} should be {}",
                x_lat,
                x_lng,
                y_lat,
                y_lng,
                output,
                expect_rect
            );
        }

        // Now we test cases where the bound does not contain nearly-antipodal
        // points, but it does contain points that are approximately 180 degrees
        // apart in latitude.
        test_rect(
            1.5,
            -FRAC_PI_2,
            1.5,
            FRAC_PI_2 - 2e-16,
            &Rect {
                lat: r1::interval::Interval::new(1.5, 1.5),
                lng: s1::interval::FULL,
            },
        );

        test_rect(
            1.5,
            -FRAC_PI_2,
            1.5,
            FRAC_PI_2 - 7e-16,
            &Rect {
                lat: r1::interval::Interval::new(1.5, 1.5),
                lng: s1::interval::Interval::new(-FRAC_PI_2, FRAC_PI_2 - 7e-16),
            },
        );
    }
}
