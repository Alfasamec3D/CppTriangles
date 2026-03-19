#include "geom.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace Geometry {

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
               const double& tolerance) {
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
               const double& tolerance) {
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
  std::cerr << "Yo, problem is not in the first if" << std::endl;

  for (int i; i < 2; ++i) n[i] = vector_product(segm_dirs[i], vect_prod_dirs);

  for (int i; i < 2; ++i)
    for (int j; j < 2; ++j)
      d[i][j] =
          approxSgn(scalar_product(n[(i + 1) % 2],
                                   segments[i]->get_points()[j] -
                                       segments[(i + 1) % 2]->get_points()[0]),
                    tolerance);
  //std::cerr << "direction of segments - "<< segm_dirs[0] << "; vector product of dirs - " <<vect_prod_dirs<<"; d - "<<d[0][0]<<std::endl;
  return ((d[0][0] != d[0][1]) && (d[1][0] != d[1][1])) ||
         ((d[0][0] == 0) && (d[0][1] == 0) &&
          ((intersect(segments[0]->get_points()[0], *segments[1], tolerance)) ||
           (intersect(segments[0]->get_points()[1], *segments[1], tolerance)) ||
           (intersect(segments[1]->get_points()[0], *segments[0], tolerance))));
}

template <typename T>
bool intersect(const VPoint<T>& vpoint, const Triangle<T>& triangle,
               const double& tolerance) {
  if (approxSgn(
          scalar_product(vpoint - triangle.get_points()[0], triangle.normal()),
          tolerance) != 0) {
    return false;
  }
  VPoint<T> triangle_normal = triangle.normal();
  VPoint<T> expected_normal;
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
               const double& tolerance) {
  VPoint<T> triangle_normal = triangle.normal();
  int d1 = approxSgn(
      scalar_product(segment.get_points()[0] - triangle.get_points()[0],
                     triangle_normal),
      tolerance);
  int d2 = approxSgn(
      scalar_product(segment.get_points()[1] - triangle.get_points()[0],
                     triangle_normal),
      tolerance);

  if ((d1 + d2 == 2) || (d1 + d2 == -2)) {
    return false;
  } else if ((d1 == 0) && (d2 == 0)) {
    bool intersect_statement = false;
    for (int i = 0; i < 3; ++i) {
      intersect_statement |= intersect(
          triangle.get_points()[(i + 1) % 3] - triangle.get_points()[i],
          segment, tolerance);
    }
    return intersect(segment.get_points()[0], triangle, tolerance) ||
           intersect(segment.get_points()[1], triangle, tolerance) ||
           intersect_statement;
  } else if ((d1 == 0) && (d2 != 0)) {
    if (intersect(segment.get_points()[0], triangle, tolerance)) {
      return true;
    }
  } else if ((d1 != 0) && (d2 == 0)) {
    if (intersect(segment.get_points()[1], triangle, tolerance)) {
      return true;
    }
  }
  return intersect(intersection(segment, triangle), triangle, tolerance);
}

template <typename T>
bool subintersect(const Triangle<T>& triangle1, const Triangle<T>& triangle2,
                  const double& tolerance) {
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
    return intersect(triangle2.get_points()[0], triangle1), tolerance;
  }
  if ((triangle1.actual_class(tolerance) == SEGMENT) &&
      (triangle2.actual_class(tolerance) == TRIANGLE)) {
    return intersect(triangle1.segment(), triangle2), tolerance;
  }

  return intersect(triangle2.segment(), triangle1), tolerance;
}

template <typename T>
bool intersect(const Triangle<T>& triangle1, const Triangle<T>& triangle2,
               const double& tolerance) {
  if ((triangle1.actual_class() == TRIANGLE) &&
      (triangle2.actual_class() == TRIANGLE)) {
    std::array<const Triangle<T>*, 2> triangle = {&triangle1, &triangle2};
    std::array<std::array<int, 3>, 2> d;
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 3; ++i) {
        d[j][i] =
            approxSgn(scalar_product(triangle[j]->get_points()[i] -
                                         triangle[(j + 1) % 2]->get_points()[0],
                                     triangle[(j + 1) % 2]->normal()),
                      tolerance);
      }
    if (((d[0][0] == d[0][1]) && (d[0][1] == d[0][2]) && (d[0][0] != 0)) ||
        ((d[1][0] == d[1][1]) && (d[1][1] == d[1][2]) && (d[1][0] != 0)))
      return false;
    if ((d[0][0] == 0) && (d[0][1] == 0) && (d[0][2] == 0))
      return intersect(triangle[0]->get_points()[0], *triangle[1], tolerance) ||
             intersect(triangle[0]->get_points()[1], *triangle[1], tolerance) ||
             intersect(triangle[0]->get_points()[2], *triangle[1], tolerance) ||
             intersect(triangle[1]->get_points()[0], *triangle[0], tolerance) ||
             intersect(triangle[1]->get_points()[1], *triangle[0], tolerance) ||
             intersect(triangle[1]->get_points()[2], *triangle[0], tolerance);
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
        // only 1 apex of triangle j lie on the plane of triangle j+1:
        if ((d[j][i] == 0) && (d[j][(i + 1) % 3] == -d[j][(i + 2) % 3])) {
          k[j] = 2;
          t[j][0] = triangle[j]->get_points()[i];
          segment[0] = {{triangle[j]->get_points()[(i + 1) % 3],
                         triangle[j]->get_points()[(i + 2) % 3]}};
          t[j][1] = intersection(segment[0], *triangle[(j + 1) % 2]);
          break;
        }

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

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

void final() {
  using numb_t = double;
  int N;
  if (!std::cin >> N) std::cerr << "invalid input for N";
  std::vector<Triangle<numb_t>> triangles;
  std::vector<AABB<numb_t>> boundingBoxes;
  std::set<int> intersectingTriangles;

  using Box = bg::model::box<bg::model::point<double, 3, bg::cs::cartesian>>;
  using Value = std::pair<Box, int>;
  bgi::rtree<Value, bgi::quadratic<16>> rtree;

  Triangle<numb_t> input_tri;
  AABB<numb_t> bound_box;
  for (int i = 0; i < N; ++i) {
    try {
      std::cin >> input_tri;
    } catch (const std::runtime_error& e) {
      std::cerr << "Problem while inputing triangle";
    }
    triangles.push_back(input_tri);

    bound_box = {input_tri.minCorner(), input_tri.maxCorner()};

    boundingBoxes.push_back(bound_box);

    Box box(bg::model::point<numb_t, 3, bg::cs::cartesian>(
                bound_box.minCorner.get_coordinates()[0],
                bound_box.minCorner.get_coordinates()[1],
                bound_box.minCorner.get_coordinates()[2]),
            bg::model::point<numb_t, 3, bg::cs::cartesian>(
                bound_box.maxCorner.get_coordinates()[0],
                bound_box.maxCorner.get_coordinates()[1],
                bound_box.maxCorner.get_coordinates()[2]));
    rtree.insert({box, i});
  }

  for (int i = 0; i < N; ++i) {
    Box queryBox(bg::model::point<numb_t, 3, bg::cs::cartesian>(
                     boundingBoxes[i].minCorner.get_coordinates()[0],
                     boundingBoxes[i].minCorner.get_coordinates()[1],
                     boundingBoxes[i].minCorner.get_coordinates()[2]),
                 bg::model::point<double, 3, bg::cs::cartesian>(
                     boundingBoxes[i].maxCorner.get_coordinates()[0],
                     boundingBoxes[i].maxCorner.get_coordinates()[1],
                     boundingBoxes[i].maxCorner.get_coordinates()[2]));

    std::vector<Value> candidates;
    rtree.query(bgi::intersects(queryBox), std::back_inserter(candidates));

    for (const auto& [box, j] : candidates) {
      if (i != j && intersect(triangles[i], triangles[j])) {
        intersectingTriangles.insert(i);
        intersectingTriangles.insert(j);
      }
    }
  }
  for (const auto& index : intersectingTriangles) {
    std::cout << index << std::endl;
  }
}

}  // namespace Geometry