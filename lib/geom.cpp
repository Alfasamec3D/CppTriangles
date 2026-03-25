#include "geom.hpp"

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace Geometry {

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

void final() {
  using numb_t = double;
  int N;
  if (!(std::cin >> N)) {
    std::cerr << "Error: Invalid input for N" << std::endl;
    return;
  }
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
      std::cerr << "Error: Problems while inputing triangle" << std::endl;
      return;
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