#include <gtest/gtest.h>

#include "geom.hpp"

TEST(ComparePointSegmentTest, BasicCase) {
  Point point{2, 1, 3};
  Point edge1{5, 5, 8};
  Point edge2{-4, -7, -7};
  Segment segment{edge1, edge2};
  bool res{true};
  EXPECT_EQ(compare(point, segment), res);
}

TEST(ComparePointTriangleTest, BasicCase) {
  Point point{2, 1, 3};
  Point apex1{5, 2, 6};
  Point apex2{-2, 2, -4};
  Point apex3{0, -2, 4};
  Triangle triangle{apex1, apex2, apex3};
  bool res{true};
  EXPECT_EQ(compare(point, triangle), res);
}

TEST(CompareSegmentsTest, BasicCase) {
  Point edge11{5, 5, 8};
  Point edge12{-4, -7, -7};
  Point edge21{14, 5, 17};
  Point edge22{-16, -5, -18};
  Segment segment1{edge11, edge12};
  Segment segment2{edge21, edge22};
  bool res{true};
  EXPECT_EQ(compare(segment1, segment2), res);
}

TEST(CompareSegmentTriangleTest, BasicCase) {
  Point edge1{5, 5, 8};
  Point edge2{-4, -7, -7};
  Point apex1{5, 2, 6};
  Point apex2{-2, 2, -4};
  Point apex3{0, -2, 4};
  Segment segment{edge1, edge2};
  Triangle triangle{apex1, apex2, apex3};
  bool res{true};
  EXPECT_EQ(compare(segment, triangle), res);
}

TEST(CompareTrianglesTest, BasicCase) {
  Point apex11{5, 2, 6};
  Point apex12{-2, 2, -4};
  Point apex13{0, -2, 4};
  Point apex21{5, 5, 8};
  Point apex22{-4, -7, -7};
  Point apex23{3, 9, -6};
  Triangle triangle1{apex11, apex12, apex13};
  Triangle triangle2{apex21, apex22, apex23};
  bool res{true};
  EXPECT_EQ(compare(triangle1, triangle2), res);
}

TEST(IntersetionTest, BasicCase) {
  Point edge1{5, 5, 8};
  Point edge2{-4, -7, -7};
  Point apex1{5, 2, 6};
  Point apex2{-2, 2, -4};
  Point apex3{0, -2, 4};
  Segment segment{edge1, edge2};
  Triangle triangle{apex1, apex2, apex3};
  Point res{2, 1, 3};
  EXPECT_EQ(intsec(segment, triangle), res);
}

TEST(ScalarProductTest, BasicCase) {
  Point point1{1, 2, 3};
  Point point2{-5, 3, -2};
  double res{-5};
  EXPECT_TRUE(eps_eq(scalprod(point1, point2), res));
}

TEST(VectorProductTest, BasicCase) {
  Vector vector1{1, 2, 3};
  Vector vector2{-5, 3, -2};
  Vector res{-13, -13, 13};
  EXPECT_EQ(vectprod(vector1, vector2), res);
}

TEST(SquaredDistanceTest, BasicCase) {
  Point point1{1, 2, 3};
  Point point2{-5, 3, -2};
  double res{62};
  EXPECT_TRUE(eps_eq(sqrdist(point1, point2), res));
}

TEST(TriangleSegmentTest, BasicCase) {
  Point point1{1, -9, 17};
  Point point2{-3, -1, 7};
  Point point3{-5, 3, 2};
  Triangle triangle{point1, point2, point3};
  Segment res{point1, point3};
  EXPECT_EQ(triangle.segm(), res);
}
