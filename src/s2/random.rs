extern crate rand;

use cgmath;
use s2::cap::Cap;
use s2::point::Point;
use s2::cellid::*;
use self::rand::Rng;

pub fn rng() -> rand::StdRng {
    rand::StdRng::new().expect("failed to get rng")
}

/// skewed_int returns a number in the range [0,2^max_log-1] with bias towards smaller numbers.
pub fn skewed_int<R>(rng: &mut R, max_log: usize) -> usize
    where R: rand::Rng
{
    let base = rng.gen_range(0, max_log + 1);
    rng.gen_range(0, 1 << 31) & ((1 << base) - 1)
}

/// cap returns a cap with a random axis such that the log of its area is
/// uniformly distributed between the logs of the two given values. The log of
/// the cap angle is also approximately uniformly distributed.
pub fn cap<R>(rng: &mut R, min_area: f64, max_area: f64) -> Cap
    where R: rand::Rng
{
    let cap_area = max_area * (min_area / max_area).powf(rng.gen_range(0., 1.));
    Cap::from_center_area(&point(rng), cap_area)
}


/// point returns a random unit-length vector.
pub fn point<R: Rng>(rng: &mut R) -> Point {
    Point::from_coords(rng.gen_range(-1., 1.),
                       rng.gen_range(-1., 1.),
                       rng.gen_range(-1., 1.))
}

pub fn frame<R: Rng>(rng: &mut R) -> cgmath::Matrix3<f64> {
    let z = point(rng);
    frame_at_point(rng, z)
}

pub fn frame_at_point<R: Rng>(rng: &mut R, z: Point) -> cgmath::Matrix3<f64> {
    let p = point(rng);
    let x = z.cross(&p).normalize();
    let y = z.cross(&x).normalize();

    cgmath::Matrix3::from_cols(x.into(), y.into(), z.into())
}

pub fn cellid<R>(rng: &mut R) -> CellID
    where R: rand::Rng
{
    let level = rng.gen_range(0, MAX_LEVEL + 1);
    cellid_for_level(rng, level)
}

pub fn cellid_for_level<R>(rng: &mut R, level: u64) -> CellID
    where R: rand::Rng
{
    let face = rng.gen_range(0, NUM_FACES as u64);
    let pos = rng.next_u64() & ((1 << POS_BITS) - 1);
    let cellid = CellID::from_face_pos_level(face, pos, level);
    assert_eq!(face, cellid.face() as u64);
    assert_eq!(level, cellid.level());

    cellid
}
