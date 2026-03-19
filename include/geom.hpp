#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <stdexcept>
#include<iostream>
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
  std::array<T, 3> coordinates_;

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
  return (-1) * vpoint;
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

/*template <typename T>
std::ostream& operator<<(std::ostream& os, const VPoint<T>& vpoint) {
  os << '(';
  for (const auto& coord : vpoint.get_coordinates()) os << coord << ' ';
  os << ')';
  return os;
}
*/
template <typename T>
class Segment {
 private:
  std::array<VPoint<T>, 2> points_;

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
  std::array<VPoint<T>, 3> points_;

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
VPoint<T> intersection(const Segment<T>& segment, const Triangle<T>& triangle);

template <typename T>
bool intersect(const VPoint<T>& vpoint, const Segment<T>& segment,
               const double& tolerance = global_tolerance);

template <typename T>
bool intersect(const Segment<T>& segment1, const Segment<T>& segment2,
               const double& tolerance = global_tolerance);

template <typename T>
bool intersect(const VPoint<T>& vpoint, const Triangle<T>& triangle,
               const double& tolerance = global_tolerance);

template <typename T>
bool intersect(const Segment<T>& segment, const Triangle<T>& triangle,
               const double& tolerance = global_tolerance);

template <typename T>
bool subintersect(const Triangle<T>& triangle1, const Triangle<T>& triangle2,
                  const double& tolerance = global_tolerance);

template <typename T>
bool intersect(const Triangle<T>& triangle1, const Triangle<T>& triangle2,
               const double& tolerance = global_tolerance);
template <typename T>
struct AABB {
  VPoint<T> minCorner, maxCorner;
};

void final();

}  // namespace Geometry