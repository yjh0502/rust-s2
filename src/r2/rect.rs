/*
Copyright 2014 Google Inc. All rights reserved.
Copyright 2017 Jihyun Yu. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

use std;

use r1::interval::{self, Interval};
use r2::point::Point;

/// Rect represents a closed axis-aligned rectangle in the (x,y) plane.
#[derive(Clone,Debug,Default)]
pub struct Rect {
    /// x interval of the rect
    pub x: Interval,
    /// y interval of the rect
    pub y: Interval,
}

/// empty rect
pub const EMPTY: Rect = Rect {
    x: interval::EMPTY,
    y: interval::EMPTY,
};

impl Rect {
    /// from_points constructs a rect that contains the given points.
    pub fn from_points(points: &[Point]) -> Self {
        // Because the default value on interval is 0,0, we need to manually
        // define the interval from the first point passed in as our starting
        // interval, otherwise we end up with the case of passing in
        // Point{0.2, 0.3} and getting the starting Rect of {0, 0.2}, {0, 0.3}
        // instead of the Rect {0.2, 0.2}, {0.3, 0.3} which is not correct.
        if points.is_empty() {
            return Self::default();
        }

        let mut r = Rect {
            x: Interval::from_point(points[0].x),
            y: Interval::from_point(points[0].y),
        };
        for p in &points[1..] {
            r = r + p;
        }
        return r;
    }

    /// from_center_size constructs a rectangle with the given center and size.
    /// Both dimensions of size must be non-negative.
    pub fn from_center_size(center: &Point, size: &Point) -> Self {
        Rect {
            x: Interval::from_point(center.x).expanded(size.x / 2.),
            y: Interval::from_point(center.y).expanded(size.y / 2.),
        }
    }

    /// empty constructs the canonical empty rectangle. Use IsEmpty() to test
    /// for empty rectangles, since they have more than one representation. A Rect{}
    /// is not the same as the EmptyRect.
    pub fn empty() -> Self {
        Rect {
            x: Interval::empty(),
            y: Interval::empty(),
        }
    }

    /// is_valid reports whether the rectangle is valid.
    /// This requires the width to be empty iff the height is empty.
    pub fn is_valid(&self) -> bool {
        self.x.is_empty() == self.y.is_empty()
    }

    /// is_empty reports whether the rectangle is empty.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// vertices returns all four vertices of the rectangle. Vertices are returned in
    /// CCW direction starting with the lower left corner.
    pub fn vertices(&self) -> [Point; 4] {
        [Point {
             x: self.x.lo,
             y: self.y.lo,
         },
         Point {
             x: self.x.hi,
             y: self.y.lo,
         },
         Point {
             x: self.x.hi,
             y: self.y.hi,
         },
         Point {
             x: self.x.lo,
             y: self.y.hi,
         }]
    }

    /// vertex_ij returns the vertex in direction i along the X-axis (0=left, 1=right) and
    /// direction j along the Y-axis (0=down, 1=up).
    pub fn vertex_ij(&self, i: isize, j: isize) -> Point {
        let x = if i == 0 { self.x.lo } else { self.x.hi };
        let y = if j == 0 { self.y.lo } else { self.y.hi };
        Point { x: x, y: y }
    }

    /// lo returns the low corner of the rect.
    pub fn lo(&self) -> Point {
        Point {
            x: self.x.lo,
            y: self.y.lo,
        }
    }

    /// hi returns the high corner of the rect.
    pub fn hi(&self) -> Point {
        Point {
            x: self.x.hi,
            y: self.y.hi,
        }
    }

    /// center returns the center of the rectangle in (x,y)-space
    pub fn center(&self) -> Point {
        Point {
            x: self.x.center(),
            y: self.y.center(),
        }
    }

    /// size returns the width and height of this rectangle in (x,y)-space. Empty
    /// rectangles have a negative width and height.
    pub fn size(&self) -> Point {
        Point {
            x: self.x.len(),
            y: self.y.len(),
        }
    }

    /// contains_point reports whether the rectangle contains the given point.
    /// Rectangles are closed regions, i.e. they contain their boundary.
    pub fn contains_point(&self, p: &Point) -> bool {
        self.x.contains(p.x) && self.y.contains(p.y)
    }

    /// interior_contains_point returns true iff the given point is contained in the interior
    /// of the region (i.e. the region excluding its boundary).
    pub fn interior_contains_point(&self, p: &Point) -> bool {
        self.x.interior_contains(p.x) && self.y.interior_contains(p.y)
    }

    /// contains reports whether the rectangle contains the given rectangle.
    pub fn contains(&self, r: &Self) -> bool {
        self.x.contains_interval(&r.x) && self.y.contains_interval(&r.y)
    }

    /// interior_contains reports whether the interior of this rectangle contains all of the
    /// points of the given other rectangle (including its boundary).
    pub fn interior_contains(&self, r: &Self) -> bool {
        self.x.interior_contains_interval(&r.x) && self.y.interior_contains_interval(&r.y)
    }

    /// intersects reports whether this rectangle and the other rectangle have any points in common.
    pub fn intersects(&self, r: &Self) -> bool {
        self.x.intersects(&r.x) && self.y.intersects(&r.y)
    }

    /// interior_intersects reports whether the interior of this rectangle intersects
    /// any point (including the boundary) of the given other rectangle.
    pub fn interior_intersects(&self, r: &Self) -> bool {
        self.x.interior_intersects(&r.x) && self.y.interior_intersects(&r.y)
    }

    /// clamp_point returns the closest point in the rectangle to the given point.
    /// The rectangle must be non-empty.
    pub fn clamp_point(&self, p: &Point) -> Point {
        Point {
            x: self.x.clamp_point(p.x),
            y: self.y.clamp_point(p.y),
        }
    }

    /// expanded returns a rectangle that has been expanded in the x-direction
    /// by margin.X, and in y-direction by margin.Y. If either margin is empty,
    /// then shrink the interval on the corresponding sides instead. The resulting
    /// rectangle may be empty. Any expansion of an empty rectangle remains empty.
    pub fn expanded(&self, margin: &Point) -> Self {
        let x = self.x.expanded(margin.x);
        let y = self.y.expanded(margin.y);
        if x.is_empty() || y.is_empty() {
            Self::empty()
        } else {
            Rect { x: x, y: y }
        }
    }

    /// expanded_by_margin returns a Rect that has been expanded by the amount on all sides.
    pub fn expanded_by_margin(&self, margin: f64) -> Self {
        self.expanded(&Point {
                          x: margin,
                          y: margin,
                      })
    }

    /// union returns the smallest rectangle containing the union of this rectangle and
    /// the given rectangle.
    pub fn union(&self, other: &Self) -> Self {
        Rect {
            x: self.x.union(&other.x),
            y: self.y.union(&other.y),
        }
    }

    /// intersection returns the smallest rectangle containing the intersection of this
    /// rectangle and the given rectangle.
    pub fn intersection(&self, other: &Self) -> Self {
        let x = self.x.intersection(&other.x);
        let y = self.y.intersection(&other.y);
        if x.is_empty() || y.is_empty() {
            Self::empty()
        } else {
            Rect { x: x, y: y }
        }
    }

    /// approx_equals returns true if the x- and y-intervals of the two rectangles are
    /// the same up to the given tolerance.
    pub fn approx_eq(&self, other: &Self) -> bool {
        self.x.approx_eq(&other.x) && self.y.approx_eq(&other.y)
    }
}

impl<'b> std::ops::Add<&'b Point> for Rect {
    type Output = Self;
    /// expands the rectangle to include the given point. The rectangle is
    /// expanded by the minimum amount possible.
    fn add(self, p: &'b Point) -> Self::Output {
        Self::Output {
            x: self.x + p.x,
            y: self.y + p.y,
        }
    }
}

impl<'b> std::ops::Add<&'b Rect> for Rect {
    type Output = Self;
    /// expands the rectangle to include the given rectangle. This is the
    /// same as replacing the rectangle by the union of the two rectangles, but
    /// is more efficient.
    fn add(self, p: &'b Rect) -> Self::Output {
        Self::Output {
            x: self.x.union(&p.x),
            y: self.y.union(&p.y),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /*
    const SW: Point = Point { x: 0., y: 0.25 };
    const SE: Point = Point { x: 0.5, y: 0.25 };
    const NE: Point = Point { x: 0.5, y: 0.75 };
    const NW: Point = Point { x: 0., y: 0.75 };
    */

    /*
    const RECT: Rect = Rect {
        x: Interval { lo: 0., hi: 0.5 },
        y: Interval {
            lo: 0.25,
            hi: 0.75,
        },
    };

    const RECT_MID: Rect = Rect {
        x: Interval {
            lo: 0.25,
            hi: 0.25,
        },
        y: Interval { lo: 0.5, hi: 0.5 },
    };

    const RECT_SW: Rect = Rect {
        x: Interval {
            lo: SW.x,
            hi: SW.x,
        },
        y: Interval {
            lo: SW.y,
            hi: SW.y,
        },
    };

    const RECT_NE: Rect = Rect {
        x: Interval {
            lo: NE.x,
            hi: NE.x,
        },
        y: Interval {
            lo: NE.y,
            hi: NE.y,
        },
    };
    */

    #[test]
    fn empty_rect() {
        assert!(EMPTY.is_valid());
        assert!(EMPTY.is_empty());
    }
}


/*
func TestFromVariousTypes(t *testing.T) {
	d1 := RectFromPoints(Point{0.1, 0}, Point{0.25, 1})
	tests := []struct {
		r1, r2 Rect
	}{
		{
			RectFromCenterSize(Point{0.3, 0.5}, Point{0.2, 0.4}),
			RectFromPoints(Point{0.2, 0.3}, Point{0.4, 0.7}),
		},
		{
			RectFromCenterSize(Point{1, 0.1}, Point{0, 2}),
			RectFromPoints(Point{1, -0.9}, Point{1, 1.1}),
		},
		{
			d1,
			Rect{d1.X, d1.Y},
		},
		{
			RectFromPoints(Point{0.15, 0.3}, Point{0.35, 0.9}),
			RectFromPoints(Point{0.15, 0.9}, Point{0.35, 0.3}),
		},
		{
			RectFromPoints(Point{0.12, 0}, Point{0.83, 0.5}),
			RectFromPoints(Point{0.83, 0}, Point{0.12, 0.5}),
		},
	}

	for _, test := range tests {
		if got := test.r1.ApproxEquals(test.r2); !got {
			t.Errorf("%v.ApproxEquals(%v); got %v want true", test.r1, test.r2, got)
		}
	}
}

func TestCenter(t *testing.T) {
	tests := []struct {
		rect Rect
		want Point
	}{
		{empty, Point{0.5, 0.5}},
		{rect, Point{0.25, 0.5}},
	}
	for _, test := range tests {
		if got := test.rect.Center(); got != test.want {
			t.Errorf("%v.Center(); got %v want %v", test.rect, got, test.want)
		}
	}
}

func TestVertices(t *testing.T) {
	want := [4]Point{sw, se, ne, nw}
	got := rect.Vertices()
	if !reflect.DeepEqual(got, want) {
		t.Errorf("%v.Vertices(); got %v want %v", rect, got, want)
	}
}

func TestContainsPoint(t *testing.T) {
	tests := []struct {
		rect Rect
		p    Point
		want bool
	}{
		{rect, Point{0.2, 0.4}, true},
		{rect, Point{0.2, 0.8}, false},
		{rect, Point{-0.1, 0.4}, false},
		{rect, Point{0.6, 0.1}, false},
		{rect, Point{rect.X.Lo, rect.Y.Lo}, true},
		{rect, Point{rect.X.Hi, rect.Y.Hi}, true},
	}
	for _, test := range tests {
		if got := test.rect.ContainsPoint(test.p); got != test.want {
			t.Errorf("%v.ContainsPoint(%v); got %v want %v", test.rect, test.p, got, test.want)
		}
	}
}

func TestInteriorContainsPoint(t *testing.T) {
	tests := []struct {
		rect Rect
		p    Point
		want bool
	}{
		// Check corners are not contained.
		{rect, sw, false},
		{rect, ne, false},
		// Check a point on the border is not contained.
		{rect, Point{0, 0.5}, false},
		{rect, Point{0.25, 0.25}, false},
		{rect, Point{0.5, 0.5}, false},
		// Check points inside are contained.
		{rect, Point{0.125, 0.6}, true},
	}
	for _, test := range tests {
		if got := test.rect.InteriorContainsPoint(test.p); got != test.want {
			t.Errorf("%v.InteriorContainsPoint(%v); got %v want %v",
				test.rect, test.p, got, test.want)
		}
	}
}

func TestIntervalOps(t *testing.T) {
	tests := []struct {
		r1, r2                                           Rect
		contains, intContains, intersects, intIntersects bool
		wantUnion, wantIntersection                      Rect
	}{
		{
			rect, rectMid,
			true, true, true, true,
			rect, rectMid,
		},
		{
			rect, rectSW,
			true, false, true, false,
			rect, rectSW,
		},
		{
			rect, rectNE,
			true, false, true, false,
			rect, rectNE,
		},
		{
			rect,
			RectFromPoints(Point{0.45, 0.1}, Point{0.75, 0.3}),
			false, false, true, true,
			RectFromPoints(Point{0, 0.1}, Point{0.75, 0.75}),
			RectFromPoints(Point{0.45, 0.25}, Point{0.5, 0.3}),
		},
		{
			rect,
			RectFromPoints(Point{0.5, 0.1}, Point{0.7, 0.3}),
			false, false, true, false,
			RectFromPoints(Point{0, 0.1}, Point{0.7, 0.75}),
			RectFromPoints(Point{0.5, 0.25}, Point{0.5, 0.3}),
		},
		{
			rect,
			RectFromPoints(Point{0.45, 0.1}, Point{0.7, 0.25}),
			false, false, true, false,
			RectFromPoints(Point{0, 0.1}, Point{0.7, 0.75}),
			RectFromPoints(Point{0.45, 0.25}, Point{0.5, 0.25}),
		},
		{
			RectFromPoints(Point{0.1, 0.2}, Point{0.1, 0.3}),
			RectFromPoints(Point{0.15, 0.7}, Point{0.2, 0.8}),
			false, false, false, false,
			RectFromPoints(Point{0.1, 0.2}, Point{0.2, 0.8}),
			EmptyRect(),
		},
		// Check that the intersection of two rectangles that overlap in x but not y
		// is valid, and vice versa.
		{
			RectFromPoints(Point{0.1, 0.2}, Point{0.4, 0.5}),
			RectFromPoints(Point{0, 0}, Point{0.2, 0.1}),
			false, false, false, false,
			RectFromPoints(Point{0, 0}, Point{0.4, 0.5}),
			EmptyRect(),
		},
		{
			RectFromPoints(Point{0, 0}, Point{0.1, 0.3}),
			RectFromPoints(Point{0.2, 0.1}, Point{0.3, 0.4}),
			false, false, false, false,
			RectFromPoints(Point{0, 0}, Point{0.3, 0.4}),
			EmptyRect(),
		},
	}
	for _, test := range tests {
		if got := test.r1.Contains(test.r2); got != test.contains {
			t.Errorf("%v.Contains(%v); got %v want %v",
				test.r1, test.r2, got, test.contains)
		}

		if got := test.r1.InteriorContains(test.r2); got != test.intContains {
			t.Errorf("%v.InteriorContains(%v); got %v want %v",
				test.r1, test.r2, got, test.contains)
		}

		if got := test.r1.Intersects(test.r2); got != test.intersects {
			t.Errorf("%v.Intersects(%v); got %v want %v",
				test.r1, test.r2, got, test.intersects)
		}

		if got := test.r1.InteriorIntersects(test.r2); got != test.intIntersects {
			t.Errorf("%v.InteriorIntersects(%v); got %v want %v",
				test.r1, test.r2, got, test.intIntersects)
		}

		tCon := test.r1.Contains(test.r2)
		if got := test.r1.Union(test.r2).ApproxEquals(test.r1); got != tCon {
			t.Errorf("%v.Union(%v) == %v.Contains(%v); got %v want %v",
				test.r1, test.r2, test.r1, test.r2, got, tCon)
		}

		tInter := test.r1.Intersects(test.r2)
		if got := !test.r1.Intersection(test.r2).IsEmpty(); got != tInter {
			t.Errorf("%v.Intersection(%v).IsEmpty() == %v.Intersects(%v); got %v want %v",
				test.r1, test.r2, test.r1, test.r2, got, tInter)
		}

		if got := test.r1.Union(test.r2); got != test.wantUnion {
			t.Errorf("%v.Union(%v); got %v want %v",
				test.r1, test.r2, got, test.wantUnion)
		}

		if got := test.r1.Intersection(test.r2); got != test.wantIntersection {
			t.Errorf("%v.Intersection(%v); got %v want %v",
				test.r1, test.r2, got, test.wantIntersection)
		}

		r := test.r1.AddRect(test.r2)

		if r != test.wantUnion {
			t.Errorf("%v.AddRect(%v); got %v want %v", test.r1, test.r2, r, test.wantUnion)
		}
	}
}

func TestAddPoint(t *testing.T) {
	r1 := rect
	r2 := EmptyRect()

	r2 = r2.AddPoint(sw)
	r2 = r2.AddPoint(se)
	r2 = r2.AddPoint(nw)
	r2 = r2.AddPoint(Point{0.1, 0.4})

	if !r1.ApproxEquals(r2) {
		t.Errorf("%v.AddPoint(%v); got false want true", r1, r2)
	}
}

func TestClampPoint(t *testing.T) {
	r := Rect{r1.Interval{Lo: 0, Hi: 0.5}, r1.Interval{Lo: 0.25, Hi: 0.75}}
	tests := []struct {
		p    Point
		want Point
	}{
		{Point{-0.01, 0.24}, Point{0, 0.25}},
		{Point{-5.0, 0.48}, Point{0, 0.48}},
		{Point{-5.0, 2.48}, Point{0, 0.75}},
		{Point{0.19, 2.48}, Point{0.19, 0.75}},

		{Point{6.19, 2.48}, Point{0.5, 0.75}},
		{Point{6.19, 0.53}, Point{0.5, 0.53}},
		{Point{6.19, -2.53}, Point{0.5, 0.25}},
		{Point{0.33, -2.53}, Point{0.33, 0.25}},
		{Point{0.33, 0.37}, Point{0.33, 0.37}},
	}
	for _, test := range tests {
		if got := r.ClampPoint(test.p); got != test.want {
			t.Errorf("%v.ClampPoint(%v); got %v want %v", r, test.p, got, test.want)
		}
	}
}

func TestExpandedEmpty(t *testing.T) {
	tests := []struct {
		rect Rect
		p    Point
	}{
		{
			EmptyRect(),
			Point{0.1, 0.3},
		},
		{
			EmptyRect(),
			Point{-0.1, -0.3},
		},
		{
			RectFromPoints(Point{0.2, 0.4}, Point{0.3, 0.7}),
			Point{-0.1, 0.3},
		},
		{
			RectFromPoints(Point{0.2, 0.4}, Point{0.3, 0.7}),
			Point{0.1, -0.2},
		},
	}
	for _, test := range tests {
		if got := test.rect.Expanded(test.p); !got.IsEmpty() {
			t.Errorf("%v.Expanded(%v); got %v want true", test.rect, test.p, got.IsEmpty())
		}
	}
}

func TestExpandedEquals(t *testing.T) {
	tests := []struct {
		rect Rect
		p    Point
		want Rect
	}{
		{
			RectFromPoints(Point{0.2, 0.4}, Point{0.3, 0.7}),
			Point{0.1, 0.3},
			RectFromPoints(Point{0.1, 0.1}, Point{0.4, 1.0}),
		},
		{
			RectFromPoints(Point{0.2, 0.4}, Point{0.3, 0.7}),
			Point{0.1, -0.1},
			RectFromPoints(Point{0.1, 0.5}, Point{0.4, 0.6}),
		},
		{
			RectFromPoints(Point{0.2, 0.4}, Point{0.3, 0.7}),
			Point{0.1, 0.1},
			RectFromPoints(Point{0.1, 0.3}, Point{0.4, 0.8}),
		},
	}
	for _, test := range tests {
		if got := test.rect.Expanded(test.p); !got.ApproxEquals(test.want) {
			t.Errorf("%v.Expanded(%v); got %v want %v", test.rect, test.p, got, test.want)
		}
	}
}
*/
