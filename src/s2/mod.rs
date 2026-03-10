mod stuv;

pub mod cell;
pub mod cellid;
pub mod cellunion;

pub mod cap;
pub mod latlng;
pub mod point;
pub mod rect;
pub mod rect_bounder;

pub mod region;

pub mod edgeutil;
pub mod metric;
pub mod predicates;

pub mod shape;

// TODO: Disable to allow testing of other modules
pub mod r#loop;

pub mod crossing_edge_query;
mod edge_clipping;
mod edge_crosser;
mod edge_crossings;
pub mod error;
pub mod lax_loop;
pub mod padded_cell;
pub mod polygon;
#[cfg(test)]
mod random;
pub mod shape_index;
pub mod shape_index_region;
mod test_util;
