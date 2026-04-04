#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace Geometry {

const double global_tolerance = 0.000001;

template <typename T>
bool approxEql(const T& object1, const T& object2,
               const double& tolerance = global_tolerance) {
  return std::abs(object1 - object2) < tolerance;
}

template <typename T>
int approxSgn(const T& object, const double& tolerance = global_tolerance) {
  if (approxEql(object, T{0}, tolerance))
    return 0;
  else if (object < 0)
    return -1;
  else
    return 1;
}

template <typename T>
class VPoint {
 private:
  std::array<T, 3> coordinates_{};

 public:
  VPoint() = default;
  VPoint(const std::array<T, 3>& coordinates) : coordinates_(coordinates) {};
  std::array<T, 3> get_coordinates() const { return {coordinates_}; }
  VPoint<T> operator+=(const VPoint<T>& vpoint) {
    for (int i = 0; i < 3; ++i) {
      coordinates_[i] += vpoint.coordinates_[i];
    }
    return *this;
  }
  template <typename U>
  VPoint<T> operator*=(const U& number) {
    for (int i = 0; i < 3; ++i) {
      coordinates_[i] *= number;
    }
    return *this;
  }

  friend std::istream& operator>>(std::istream& is, VPoint<T>& vpoint) {
    for (auto& coord : vpoint.coordinates_)
      if (!(is >> coord))
        throw std::runtime_error(
            "Fatal: input was not valid for coordinate of point");
    return is;
  }
};

template <typename T>
VPoint<T> operator+(const VPoint<T>& vpoint1, const VPoint<T>& vpoint2) {
  VPoint<T> copy_vpoint1 = vpoint1;
  return copy_vpoint1 += vpoint2;
}

template <typename T, typename U>
VPoint<T> operator*(const U& number, const VPoint<T>& vpoint) {
  VPoint<T> copy_vpoint = vpoint;
  return copy_vpoint *= number;
}

template <typename T, typename U>
VPoint<T> operator*(const VPoint<T>& vpoint, const U& number) {
  return number * vpoint;
}

template <typename T, typename U>
VPoint<T> operator/(const VPoint<T>& vpoint, const U& number) {
  return vpoint * (1 / number);
}

template <typename T>
VPoint<T> operator-(const VPoint<T>& vpoint) {
  return T{-1} * vpoint;
}

template <typename T>
VPoint<T> operator-(const VPoint<T>& vpoint1, const VPoint<T>& vpoint2) {
  return vpoint1 + -vpoint2;
}

template <typename T>
bool approxEql(const VPoint<T>& object1, const VPoint<T>& object2,
               double tolerance = global_tolerance) {
  VPoint<T> object_difference = object1 - object2;
  for (const auto& coordinate : object_difference.get_coordinates())
    if (approxSgn(coordinate, tolerance) != 0) return false;
  return true;
}

template <typename T>
T scalar_product(const VPoint<T>& point1, const VPoint<T>& point2) {
  T result{};
  for (int i = 0; i < 3; ++i)
    result += point1.get_coordinates()[i] * point2.get_coordinates()[i];

  return result;
}

template <typename T>
VPoint<T> vector_product(const VPoint<T>& point1, const VPoint<T>& point2) {
  std::array<T, 3> coordinates{};
  for (int i = 0; i < 3; ++i)
    coordinates[i] = point1.get_coordinates()[(i + 1) % 3] *
                         point2.get_coordinates()[(i + 2) % 3] -
                     point1.get_coordinates()[(i + 2) % 3] *
                         point2.get_coordinates()[(i + 1) % 3];
  return {coordinates};
}

template <typename T>
T sqr_distance(const VPoint<T>& point1, const VPoint<T>& point2) {
  const VPoint<T> difference_point = point1 - point2;

  return scalar_product(difference_point, difference_point);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const VPoint<T>& vpoint) {
  os << '(';
  for (const auto& coord : vpoint.get_coordinates()) os << coord << ' ';
  os << ')';
  return os;
}

template <typename T>
class Segment {
 private:
  std::array<VPoint<T>, 2> points_{};

 public:
  Segment() = default;
  Segment(const std::array<VPoint<T>, 2>& points) : points_(points) {};

  std::array<VPoint<T>, 2> get_points() const { return points_; }

  VPoint<T> direction() const { return points_[1] - points_[0]; }
};

template <typename T>
bool approxEql(const Segment<T>& object1, const Segment<T>& object2,
               double tolerance = global_tolerance) {
  return approxEql(object1.get_points()[0], object2.get_points()[0],
                   tolerance) &&
             approxEql(object1.get_points()[1], object2.get_points()[1],
                       tolerance) ||
         approxEql(object1.get_points()[0], object2.get_points()[1],
                   tolerance) &&
             approxEql(object1.get_points()[1], object2.get_points()[0],
                       tolerance);
}

enum GeometryClass { VPOINT, SEGMENT, TRIANGLE };
template <typename T>
class Triangle {
 private:
  std::array<VPoint<T>, 3> points_{};

 public:
  Triangle() = default;
  Triangle(const std::array<VPoint<T>, 3>& points) : points_(points) {};

  std::array<VPoint<T>, 3> get_points() const { return points_; }

  VPoint<T> normal() const {
    return vector_product(points_[1] - points_[0], points_[2] - points_[1]);
  }

  GeometryClass actual_class(const double& tolerance = global_tolerance) const {
    if (approxEql(points_[0], points_[1], tolerance) &&
        approxEql(points_[1], points_[2], tolerance))
      return VPOINT;
    if (approxEql(normal(), {{0, 0, 0}})) return SEGMENT;
    return TRIANGLE;
  }

  Segment<T> segment() const {
    assert(actual_class() != TRIANGLE);
    int k;
    for (int i = 0; i < 3; ++i) {
      if ((sqr_distance(points_[i], points_[(i + 1) % 3]) >=
           sqr_distance(points_[(i + 1) % 3], points_[(i + 2) % 3])) &&
          (sqr_distance(points_[i], points_[(i + 1) % 3]) >=
           sqr_distance(points_[(i + 2) % 3], points_[i]))) {
        k = i;
        break;
      }
    }
    return {{points_[k], points_[(k + 1) % 3]}};
  }

  friend std::istream& operator>>(std::istream& is, Triangle<T>& triangle) {
    for (auto& point : triangle.points_) is >> point;
    return is;
  }

  VPoint<T> minCorner() const {
    std::array<T, 3> coordinates;
    for (int i = 0; i < 3; ++i)
      coordinates[i] = std::min({points_[0].get_coordinates()[i],
                                 points_[1].get_coordinates()[i],
                                 points_[2].get_coordinates()[i]});
    return {coordinates};
  }
  VPoint<T> maxCorner() const {
    std::array<T, 3> coordinates;
    for (int i = 0; i < 3; ++i)
      coordinates[i] = std::max({points_[0].get_coordinates()[i],
                                 points_[1].get_coordinates()[i],
                                 points_[2].get_coordinates()[i]});
    return {coordinates};
  }
};

template <typename T>
VPoint<T> intersection(const Segment<T>& segment, const Triangle<T>& triangle) {
  return segment.get_points()[0] +
         segment.direction() *
             scalar_product(triangle.get_points()[0] - segment.get_points()[0],
                            triangle.normal()) /
             scalar_product(segment.direction(), triangle.normal());
}

template <typename T>
bool intersect(const VPoint<T>& vpoint, const Segment<T>& segment,
               const double& tolerance = global_tolerance) {
  auto scal_prod = scalar_product(vpoint - segment.get_points()[0],
                                  vpoint - segment.get_points()[1]);
  return approxEql(std::pow(scal_prod, 2),
                   sqr_distance(segment.get_points()[0], vpoint) *
                       sqr_distance(segment.get_points()[1], vpoint),
                   tolerance) &&
         (approxSgn(scal_prod, tolerance) <= 0);
}

template <typename T>
bool intersect(const Segment<T>& segment1, const Segment<T>& segment2,
               const double& tolerance = global_tolerance) {
  const std::array<const Segment<T>*, 2> segments = {&segment1, &segment2};
  const std::array<VPoint<T>, 2> segm_dirs = {segments[0]->direction(),
                                              segments[1]->direction()};
  const VPoint<T> vect_prod_dirs = vector_product(segm_dirs[0], segm_dirs[1]);

  std::array<VPoint<T>, 2> n;
  std::array<std::array<int, 2>, 2> d;

  if (approxSgn(scalar_product(
                    segments[1]->get_points()[0] - segments[0]->get_points()[0],
                    vect_prod_dirs),
                tolerance) != 0)
    return false;

  for (int i = 0; i < 2; ++i)
    n[i] = vector_product(segm_dirs[i], vect_prod_dirs);

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      d[i][j] =
          approxSgn(scalar_product(n[(i + 1) % 2],
                                   segments[i]->get_points()[j] -
                                       segments[(i + 1) % 2]->get_points()[0]),
                    tolerance);

  return ((d[0][0] != d[0][1]) && (d[1][0] != d[1][1])) ||
         ((d[0][0] == 0) && (d[0][1] == 0) &&
          ((intersect(segments[0]->get_points()[0], *segments[1], tolerance)) ||
           (intersect(segments[0]->get_points()[1], *segments[1], tolerance)) ||
           (intersect(segments[1]->get_points()[0], *segments[0], tolerance))));
}

template <typename T>
bool intersect(const VPoint<T>& vpoint, const Triangle<T>& triangle,
               const double& tolerance = global_tolerance) {
  if (approxSgn(
          scalar_product(vpoint - triangle.get_points()[0], triangle.normal()),
          tolerance) != 0) {
    return false;
  }
  VPoint<T> triangle_normal = triangle.normal();
  VPoint<T> expected_normal{};
  VPoint<T> custom_normal;
  Triangle<T> custom_triangle;
  for (int i = 0; i < 3; ++i) {
    custom_triangle = Triangle<T>{
        {vpoint, triangle.get_points()[i], triangle.get_points()[(i + 1) % 3]}};
    custom_normal = custom_triangle.normal();
    expected_normal +=
        custom_normal *
        approxSgn(scalar_product(triangle_normal, custom_normal), tolerance);
  }
  return approxEql(triangle_normal, expected_normal, tolerance);
}

template <typename T>
bool intersect(const Segment<T>& segment, const Triangle<T>& triangle,
               const double& tolerance = global_tolerance) {
  VPoint<T> triangle_normal = triangle.normal();
  int d1 = approxSgn(
      scalar_product(segment.get_points()[0] - triangle.get_points()[0],
                     triangle_normal),
      tolerance);
  int d2 = approxSgn(
      scalar_product(segment.get_points()[1] - triangle.get_points()[0],
                     triangle_normal),
      tolerance);

  // segment is above triangle:
  if ((d1 + d2 == 2) || (d1 + d2 == -2)) {
    return false;
    // segment is on plane of triangle:
  } else if ((d1 == 0) && (d2 == 0))
    return intersect(segment.get_points()[0], triangle, tolerance) ||
           intersect(segment.get_points()[1], triangle, tolerance) ||
           intersect(
               segment,
               Segment<T>{{triangle.get_points()[0], triangle.get_points()[1]}},
               tolerance) ||
           intersect(
               segment,
               Segment<T>{{triangle.get_points()[1], triangle.get_points()[2]}},
               tolerance) ||
           intersect(
               segment,
               Segment<T>{{triangle.get_points()[2], triangle.get_points()[0]}},
               tolerance);
  // only the first point of segment is on plane of triangle:
  else if ((d1 == 0) && (d2 != 0))
    return intersect(segment.get_points()[0], triangle, tolerance);
  // only the second point of segment is on plane of triangle:
  else if ((d1 != 0) && (d2 == 0))
    return intersect(segment.get_points()[1], triangle, tolerance);
  // points of segment are seperated by plane of triangles:
  return intersect(intersection(segment, triangle), triangle, tolerance);
}

template <typename T>
bool subintersect(const Triangle<T>& triangle1, const Triangle<T>& triangle2,
                  const double& tolerance = global_tolerance) {
  if ((triangle1.actual_class(tolerance) == VPOINT) &&
      (triangle2.actual_class(tolerance) == VPOINT)) {
    return approxEql(triangle1.get_points()[0], triangle2.get_points()[0],
                     tolerance);
  }
  if ((triangle1.actual_class(tolerance) == VPOINT) &&
      (triangle2.actual_class(tolerance) == SEGMENT)) {
    return intersect(triangle1.get_points()[0], triangle2.segment(), tolerance);
  }
  if ((triangle1.actual_class(tolerance) == SEGMENT) &&
      (triangle2.actual_class(tolerance) == VPOINT)) {
    return intersect(triangle2.get_points()[0], triangle1.segment(), tolerance);
  }
  if ((triangle1.actual_class(tolerance) == SEGMENT) &&
      (triangle2.actual_class(tolerance) == SEGMENT)) {
    return intersect(triangle1.segment(), triangle2.segment(), tolerance);
  }
  if ((triangle1.actual_class(tolerance) == VPOINT) &&
      (triangle2.actual_class(tolerance) == TRIANGLE)) {
    return intersect(triangle1.get_points()[0], triangle2, tolerance);
  }
  if ((triangle1.actual_class(tolerance) == TRIANGLE) &&
      (triangle2.actual_class(tolerance) == VPOINT)) {
    return intersect(triangle2.get_points()[0], triangle1, tolerance);
  }
  if ((triangle1.actual_class(tolerance) == SEGMENT) &&
      (triangle2.actual_class(tolerance) == TRIANGLE)) {
    return intersect(triangle1.segment(), triangle2, tolerance);
  }

  return intersect(triangle2.segment(), triangle1, tolerance);
}

template <typename T>
bool intersect(const Triangle<T>& triangle1, const Triangle<T>& triangle2,
               const double& tolerance = global_tolerance) {
  if ((triangle1.actual_class() == TRIANGLE) &&
      (triangle2.actual_class() == TRIANGLE)) {
    std::array<const Triangle<T>*, 2> triangle = {&triangle1, &triangle2};
    std::array<std::array<int, 3>, 2> d;

    // set values for d:
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 3; ++i) {
        d[j][i] =
            approxSgn(scalar_product(triangle[j]->get_points()[i] -
                                         triangle[(j + 1) % 2]->get_points()[0],
                                     triangle[(j + 1) % 2]->normal()),
                      tolerance);
      }

    // both triangles lie above each other:
    if (((d[0][0] == d[0][1]) && (d[0][1] == d[0][2]) && (d[0][0] != 0)) ||
        ((d[1][0] == d[1][1]) && (d[1][1] == d[1][2]) && (d[1][0] != 0)))
      return false;

    // both triangles are on the same plane:
    if ((d[0][0] == 0) && (d[0][1] == 0) && (d[0][2] == 0))
      return intersect(triangle[0]->get_points()[0], *triangle[1], tolerance) ||
             intersect(triangle[0]->get_points()[1], *triangle[1], tolerance) ||
             intersect(triangle[0]->get_points()[2], *triangle[1], tolerance) ||
             intersect(triangle[1]->get_points()[0], *triangle[0], tolerance) ||
             intersect(triangle[1]->get_points()[1], *triangle[0], tolerance) ||
             intersect(triangle[1]->get_points()[2], *triangle[0], tolerance) ||
             intersect(Segment<T>{{triangle[0]->get_points()[0],
                                   triangle[0]->get_points()[1]}},
                       *triangle[1], tolerance);

    // Having examined these 2 cases, we understand, the polygons intersect each
    // other's planes, but do not lie on them
    std::array<Segment<T>, 2> segment;
    std::array<std::array<VPoint<T>, 2>, 2> t;
    std::array<int, 2> k;
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 3; ++i) {
        // every apex of triangle j doesn't lie on the triangle j+1's plane:
        if ((d[j][i] == -d[j][(i + 1) % 3]) &&
            (d[j][(i + 1) % 3] == d[j][(i + 2) % 3])) {
          k[j] = 2;
          segment[0] = {{triangle[j]->get_points()[i],
                         triangle[j]->get_points()[(i + 1) % 3]}};
          segment[1] = {{triangle[j]->get_points()[(i + 2) % 3],
                         triangle[j]->get_points()[i]}};
          t[j][0] = intersection(segment[0], *triangle[(j + 1) % 2]);
          t[j][1] = intersection(segment[1], *triangle[(j + 1) % 2]);
          break;
        }
        // only 2 apexes of triangle j lie on the plane of triangle j+1:
        if ((d[j][i] == 0) && (d[j][(i + 1) % 3] == 0)) {
          k[j] = 2;
          t[j][0] = triangle[j]->get_points()[i];
          t[j][1] = triangle[j]->get_points()[(i + 1) % 3];
          break;
        }
        // 1 apex of triangle j lies on the plane of triangle j+1, 2 others lie
        // on opposite sides of space, separated by plane of triangle j+1:
        if ((d[j][i] == 0) && (d[j][(i + 1) % 3] == -d[j][(i + 2) % 3])) {
          k[j] = 2;
          t[j][0] = triangle[j]->get_points()[i];
          segment[0] = {{triangle[j]->get_points()[(i + 1) % 3],
                         triangle[j]->get_points()[(i + 2) % 3]}};
          t[j][1] = intersection(segment[0], *triangle[(j + 1) % 2]);
          break;
        }

        // 1 apex of triangle j lies on the plane of triangle j+1, 2 others lie
        // on the same side of space, separated by plane of triangle j+1:
        if ((d[j][i] == 0) && (d[j][(i + 1) % 3] == d[j][(i + 2) % 3])) {
          k[j] = 1;
          t[j][0] = triangle[j]->get_points()[i];
          break;
        }
      }
    std::array<Segment<T>, 2> final_segments = {Segment<T>{{t[1][0], t[1][1]}},
                                                Segment<T>{{t[0][0], t[0][1]}}};
    if ((k[0] == 1) && (k[1] == 1)) {
      return approxEql(t[0][1], t[0][0], tolerance);
    } else if ((k[0] == 1) && (k[1] == 2)) {
      return intersect(t[0][0], final_segments[0], tolerance);
    } else if ((k[0] == 2) && (k[1] == 1)) {
      return intersect(t[1][0], final_segments[1], tolerance);
    } else {
      return intersect(final_segments[0], final_segments[1], tolerance);
    }
  }
  // if one of the triangles is actually not a triangle
  return subintersect(triangle1, triangle2, tolerance);
}

template <typename T>
struct AABB {
  VPoint<T> minCorner, maxCorner;
};

void final();

}  // namespace Geometry