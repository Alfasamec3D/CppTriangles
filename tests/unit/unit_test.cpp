#include <gtest/gtest.h>

#include "geom.hpp"
using namespace Geometry;

using numb_t = double;

TEST(IntersectPointSegmentTest, BasicCase) {
  VPoint<numb_t> point{{2, 1, 3}};
  VPoint<numb_t> edge1{{5, 5, 8}};
  VPoint<numb_t> edge2{{-4, -7, -7}};
  Segment<numb_t> segment{{edge1, edge2}};
  bool res{true};
  EXPECT_EQ(intersect(point, segment), res);
}

TEST(IntersectPointTriangleTest, BasicCase) {
  VPoint<numb_t> point{{2, 1, 3}};
  VPoint<numb_t> apex1{{5, 2, 6}};
  VPoint<numb_t> apex2{{-2, 2, -4}};
  VPoint<numb_t> apex3{{0, -2, 4}};
  Triangle<numb_t> triangle{{apex1, apex2, apex3}};
  bool res{true};
  EXPECT_EQ(intersect(point, triangle), res);
}

TEST(IntersectSegmentsTest, BasicCase) {
  VPoint<numb_t> edge11{{5, 5, 8}};
  VPoint<numb_t> edge12{{-4, -7, -7}};
  VPoint<numb_t> edge21{{14, 5, 17}};
  VPoint<numb_t> edge22{{-16, -5, -18}};
  Segment<numb_t> segment1{{edge11, edge12}};
  Segment<numb_t> segment2{{edge21, edge22}};
  bool res{true};
  EXPECT_EQ(intersect(segment1, segment2), res);
}

TEST(IntersectSegmentTriangleTest, BasicCase) {
  VPoint<numb_t> edge1{{5, 5, 8}};
  VPoint<numb_t> edge2{{-4, -7, -7}};
  VPoint<numb_t> apex1{{5, 2, 6}};
  VPoint<numb_t> apex2{{-2, 2, -4}};
  VPoint<numb_t> apex3{{0, -2, 4}};
  Segment<numb_t> segment{{edge1, edge2}};
  Triangle<numb_t> triangle{{apex1, apex2, apex3}};
  bool res{true};
  EXPECT_EQ(intersect(segment, triangle), res);
}

TEST(IntersectTrianglesTest, BasicCase) {
  VPoint<numb_t> apex11{{5, 2, 6}};
  VPoint<numb_t> apex12{{-2, 2, -4}};
  VPoint<numb_t> apex13{{0, -2, 4}};
  VPoint<numb_t> apex21{{5, 5, 8}};
  VPoint<numb_t> apex22{{-4, -7, -7}};
  VPoint<numb_t> apex23{{3, 9, -6}};
  Triangle<numb_t> triangle1{{apex11, apex12, apex13}};
  Triangle<numb_t> triangle2{{apex21, apex22, apex23}};
  bool res{true};
  EXPECT_EQ(intersect(triangle1, triangle2), res);
}

TEST(IntersetionTest, BasicCase) {
  VPoint<numb_t> edge1{{5, 5, 8}};
  VPoint<numb_t> edge2{{-4, -7, -7}};
  VPoint<numb_t> apex1{{5, 2, 6}};
  VPoint<numb_t> apex2{{-2, 2, -4}};
  VPoint<numb_t> apex3{{0, -2, 4}};
  Segment<numb_t> segment{{edge1, edge2}};
  Triangle<numb_t> triangle{{apex1, apex2, apex3}};
  VPoint<numb_t> res{{2, 1, 3}};
  EXPECT_TRUE(approxEql(intersection(segment, triangle), res));
}

TEST(ScalarProductTest, BasicCase) {
  VPoint<numb_t> point1{{1, 2, 3}};
  VPoint<numb_t> point2{{-5, 3, -2}};
  double res{-5};
  EXPECT_TRUE(approxEql(scalar_product(point1, point2), res));
}

TEST(VectorProductTest, BasicCase) {
  VPoint<numb_t> vector1{{1, 2, 3}};
  VPoint<numb_t> vector2{{-5, 3, -2}};
  VPoint<numb_t> res{{-13, -13, 13}};
  EXPECT_TRUE(approxEql(vector_product(vector1, vector2), res));
}

TEST(SquaredDistanceTest, BasicCase) {
  VPoint<numb_t> point1{{1, 2, 3}};
  VPoint<numb_t> point2{{-5, 3, -2}};
  double res{62};
  EXPECT_TRUE(approxEql(sqr_distance(point1, point2), res));
}

TEST(TriangleSegmentTest, BasicCase) {
  VPoint<numb_t> point1{{1, -9, 17}};
  VPoint<numb_t> point2{{-3, -1, 7}};
  VPoint<numb_t> point3{{-5, 3, 2}};
  Triangle<numb_t> triangle{{point1, point2, point3}};
  Segment<numb_t> res{{point1, point3}};
  EXPECT_TRUE(approxEql(triangle.segment(), res));
}

TEST(DirectionTest, BasicCase) {
  VPoint<numb_t> edge1{{5, 5, 8}};
  VPoint<numb_t> edge2{{-4, -7, -7}};
  Segment<numb_t> segment{{edge1, edge2}};
  VPoint<numb_t> expected_direction{{-9, -12, -15}};
  EXPECT_TRUE(approxEql(segment.direction(), expected_direction));
}

TEST(ApproxSgnTest, BasicCase) {
  numb_t value = -0.001;
  int expected_sgn = -1;
  int real_sgn = approxSgn(value, 0.0001);
  EXPECT_EQ(real_sgn, expected_sgn);
};

TEST(VPointDifferenceTest, BasicCase) {
  VPoint<numb_t> vpoint1{{5, 5, 8}};
  VPoint<numb_t> vpoint2{{-4, -7, -7}};
  VPoint<numb_t> real_difference = vpoint2 - vpoint1;
  VPoint<numb_t> expected_difference{{-9, -12, -15}};
  EXPECT_TRUE(approxEql(real_difference, expected_difference));
}

TEST(GetPointsTest, BasicCase) {
  VPoint<numb_t> vpoint1{{5, 5, 8}};
  VPoint<numb_t> vpoint2{{-4, -7, -7}};
  Segment<numb_t> segment{{vpoint1, vpoint2}};
  std::array<VPoint<numb_t>, 2> real_points = segment.get_points();
  std::array<VPoint<numb_t>, 2> expected_points = {vpoint1, vpoint2};
  for (int i = 0; i < 2; ++i)
    EXPECT_TRUE(approxEql(real_points[i], expected_points[i]));
}

TEST(ApproxEqlNumberTest, BasicCase) {
  numb_t a = 2.000001;
  numb_t b = 2.000001;
  bool result = approxEql(a, b, 0.0000001);
  EXPECT_TRUE(result);
}

TEST(ApproxEqlVpointTest, BasicCase) {
  VPoint<numb_t> vpoint1{{5.0002, 5.00003, 8.004}};
  VPoint<numb_t> vpoint2{{5.0002, 5.00003, 8.004}};
  bool result = approxEql(vpoint1, vpoint2, 0.000001);
  EXPECT_TRUE(result);
}

TEST(ApproxEqlSegmentTest, BasicCase) {
  VPoint<numb_t> vpoint1{{5.0002, 5.00003, 8.004}};
  VPoint<numb_t> vpoint2{{-4.001, -7.0006, -7.0008}};
  Segment<numb_t> segment1{{vpoint1, vpoint2}};
  Segment<numb_t> segment2{{vpoint2, vpoint1}};
  EXPECT_TRUE(approxEql(segment1, segment2, 0.000001));
}