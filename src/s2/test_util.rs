//
// // random_uniform_int returns a uniformly distributed integer in the range [lo,hi).
// pub fn random_uniform_int(lo: isize, hi: isize) -> isize {
//     if lo >= hi {
//         return lo;
//     }
//     lo + (rand::random::<isize>() % (hi - lo))
// }
//
// // Legacy function that returns a uniformly distributed integer in the range [0,n).
// pub fn random_uniform_int_legacy(bound: isize) -> isize {
//     random::<isize>() * (bound as f64)
// }
// pub fn random_uniform_usize() -> usize {
//     rand::random::<usize>()
// }
//
// // random_uniform_float64 returns a uniformly distributed value in the range [min, max).
// pub fn random_uniform_float64(min: f64, max: f64) -> f64 {
//     min + rand::random::<f64>() * (max - min)
// }
//
// // random_float64 returns a uniformly distributed value in the range [0, 1).
// pub fn random_float64() -> f64 {
//     rand::random::<f64>()
// }
//
// // one_in returns true with a probability of 1/n.
// pub fn one_in(n: isize) -> bool {
//     rand::random::<isize>() % n == 0
// }
//
// // one_in_n returns true with a probability of 1/n.
// pub fn one_in_n(n: u32) -> bool {
//     rand::random::<u32>() % n == 0
// }
//
// // randomCellID returns a random CellID at a randomly chosen
// // level. The distribution is uniform over the space of cell ids,
// // but only approximately uniform over the surface of the sphere.
// pub fn random_cell_id() -> CellID {
//     random_cell_id_for_level(random_uniform_int_legacy((MAX_LEVEL + 1) as isize))
// }
//
// // randomCellIDForLevel returns a random CellID at the given level.
// // The distribution is uniform over the space of cell ids, but only
// // approximately uniform over the surface of the sphere.
// pub fn random_cell_id_for_level(level: isize) -> CellID {
//     let face = random_uniform_int_legacy(NUM_FACES as isize);
//     let pos = random_uniform_usize() & ((1_isize.wrapping_shl(POS_BITS as u32)) - 1) as usize;
//     CellID::from_face_pos_level(face as u64, pos as u64, level as u64)
// }

#[macro_export]
macro_rules! assert_eq_with_callout {
    // Basic version with default assertion text
    ($x:expr, $y:expr) => {
        assert_eq_with_callout!($x, $y, "assertion failed: ");
    };

    // With customerror message
    ($x:expr, $y:expr, $msg:expr) => {
        let x_val = $x;
        let y_val = $y;
        let diff = x_val != y_val;
        // let difference = x_val - y_val;
        if diff {
            panic!(
                concat!(
                    "{}\n",
                    "Left: {:?}\n",
                    "Right:   {:?}\n",
                    "NOT EQUAL!",
                    // "Diff:     {}\n",
                ),
                $msg,
                x_val,
                y_val,
                // difference
            );
        }
    };
}
