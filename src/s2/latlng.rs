
use std;
use std::f64::consts::PI;
use consts::remainder;
use r3::vector::Vector;
use s1::angle::*;
use s2::point::Point;

const NORTH_POLE_LAT: Angle = Angle(PI / 2.);
const SOUTH_POLE_LAT: Angle = Angle(PI / -2.);

#[derive(Clone)]
pub struct LatLng {
    pub lat: Angle,
    pub lng: Angle,
}

impl std::fmt::Debug for LatLng {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[{:?}, {:?}]", self.lat, self.lng)
    }
}

impl LatLng {
    pub fn from_degrees(lat: Deg, lng: Deg) -> Self {
        LatLng {
            lat: lat.into(),
            lng: lng.into(),
        }
    }

    pub fn is_valid(&self) -> bool {
        self.lat.0.abs() <= PI / 2. && self.lng.0.abs() <= PI
    }

    pub fn normalized(&self) -> Self {
        let lat = if self.lat.0 > NORTH_POLE_LAT.0 {
            NORTH_POLE_LAT
        } else if self.lat.0 < SOUTH_POLE_LAT.0 {
            SOUTH_POLE_LAT
        } else {
            self.lat.clone()
        };

        LatLng {
            lat: lat,
            lng: Angle(remainder(self.lng.0, PI * 2.)),
        }
    }

    pub fn distance(&self, other: &Self) -> Angle {
        let dlat = (0.5 * (other.lat.0 - self.lat.0)).sin();
        let dlng = (0.5 * (other.lng.0 - self.lng.0)).sin();

        let x = dlat * dlat + dlng * dlng * self.lat.0.cos() * other.lat.0.cos();
        Angle(2. * x.sqrt().atan2((1. - x).max(0.).sqrt()))
    }
}

impl Point {
    pub fn latitude(&self) -> Angle {
        let v = &self.0;
        let l = (v.x * v.x + v.y * v.y).sqrt();
        Angle(v.z.atan2(l))
    }

    pub fn longitude(&self) -> Angle {
        let v = &self.0;
        Angle(v.y.atan2(v.x))
    }
}

impl From<LatLng> for Point {
    fn from(ll: LatLng) -> Self {
        let phi = ll.lat.0;
        let theta = ll.lng.0;
        let cosphi = phi.cos();
        Point(Vector {
                  x: theta.cos() * cosphi,
                  y: theta.sin() * cosphi,
                  z: phi.sin(),
              })
    }
}

impl<'a> From<&'a Point> for LatLng {
    fn from(p: &'a Point) -> Self {
        LatLng {
            lat: p.latitude(),
            lng: p.longitude(),
        }
    }
}
impl From<Point> for LatLng {
    fn from(p: Point) -> Self {
        LatLng::from(&p)
    }
}

#[cfg(test)]
mod tests {
    use s1;
    use r3::vector::Vector;
    use s2::point::Point;
    use super::*;

    macro_rules! ll {
        ($lat: expr, $lng: expr) => {
            LatLng{lat: s1::angle::Deg($lat).into(), lng: s1::angle::Deg($lng).into()}
        }
    }
    macro_rules! p {
        ($x: expr, $y: expr, $z: expr) => {
            Point(Vector{x:$x as f64, y:$y as f64, z:$z as f64})
        }
    }

    fn test_latlng_normalized_case(descr: &str, pos: LatLng, want: LatLng) {
        let desc: String = descr.into();
        let normalized = pos.normalized();
        assert!(normalized.is_valid(), desc);

        let distance = normalized.distance(&want);
        assert!(distance < s1::angle::Deg(1e-13).into(), desc);
    }

    #[test]
    fn test_latlng_normalized() {
        test_latlng_normalized_case(&"Valid lat/lng",
                                    ll!(21.8275043, 151.1979675),
                                    ll!(21.8275043, 151.1979675));
        test_latlng_normalized_case(&"Valid lat/lng in the West",
                                    ll!(21.8275043, -151.1979675),
                                    ll!(21.8275043, -151.1979675));
        test_latlng_normalized_case(&"Beyond the North pole",
                                    ll!(95., 151.1979675),
                                    ll!(90., 151.1979675));
        test_latlng_normalized_case("Beyond the South pole",
                                    ll!(-95., 151.1979675),
                                    ll!(-90., 151.1979675));
        test_latlng_normalized_case("At the date line (from East)",
                                    ll!(21.8275043, 180.),
                                    ll!(21.8275043, 180.));
        test_latlng_normalized_case("At the date line (from West)",
                                    ll!(21.8275043, -180.),
                                    ll!(21.8275043, -180.));
        test_latlng_normalized_case("Across the date line going East",
                                    ll!(21.8275043, 181.0012),
                                    ll!(21.8275043, -178.9988));
        test_latlng_normalized_case("Across the date line going West",
                                    ll!(21.8275043, -181.0012),
                                    ll!(21.8275043, 178.9988));
        test_latlng_normalized_case("All wrong", ll!(256., 256.), ll!(90., -104.));
    }

    fn test_approx_eq(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-14);
    }

    fn test_latlng_point_conversion_case(ll: LatLng, p: Point) {
        //TODO
        let llp: Point = ll.clone().into();
        test_approx_eq(llp.0.x, p.0.x);
        test_approx_eq(llp.0.y, p.0.y);
        test_approx_eq(llp.0.z, p.0.z);

        let pll: LatLng = p.into();
        test_approx_eq(pll.lat.0, ll.lat.0);
        let is_polar = ll.lng.0 == PI / 2. || ll.lng.0 == PI / -2.;
        if !is_polar {
            test_approx_eq(pll.lng.0, ll.lng.0);
        }
    }

    #[test]
    fn test_latlng_point_conversion() {
        test_latlng_point_conversion_case(ll!(0., 0.), p!(1, 0, 0));
        test_latlng_point_conversion_case(ll!(90., 0.), p!(6.12323e-17, 0, 1));
        test_latlng_point_conversion_case(ll!(-90., 0.), p!(6.12323e-17, 0, -1));
        test_latlng_point_conversion_case(ll!(0., 180.), p!(-1, 1.22465e-16, 0));
        test_latlng_point_conversion_case(ll!(0., -180.), p!(-1, -1.22465e-16, 0));
        test_latlng_point_conversion_case(ll!(90., 180.), p!(-6.12323e-17, 7.4988e-33, 1));
        test_latlng_point_conversion_case(ll!(90., -180.), p!(-6.12323e-17, -7.4988e-33, 1));
        test_latlng_point_conversion_case(ll!(-90., 180.), p!(-6.12323e-17, 7.4988e-33, -1));
        test_latlng_point_conversion_case(ll!(-90., -180.), p!(-6.12323e-17, -7.4988e-33, -1));
        test_latlng_point_conversion_case(ll!(-81.82750430354997, 151.19796752929685),
                                          p!(-0.12456788151479525,
                                             0.0684875268284729,
                                             -0.989844584550441));
    }

    fn test_latlng_distance_case(ll1: LatLng, ll2: LatLng, want: f64, tolerance: f64) {
        let distance: s1::angle::Deg = ll1.distance(&ll2).into();
        assert!((distance.0 - want).abs() <= tolerance);
    }

    #[test]
    fn test_latlng_distance() {
        test_latlng_distance_case(ll!(90., 0.), ll!(90., 0.), 0., 0.);
        test_latlng_distance_case(ll!(-37., 25.), ll!(-66., -155.), 77., 1e-13);
        test_latlng_distance_case(ll!(0., 165.), ll!(0., -80.), 115., 1e-13);
        test_latlng_distance_case(ll!(47., -127.), ll!(-47., 53.), 180., 2e-6);
    }
}
