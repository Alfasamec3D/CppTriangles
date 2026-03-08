#include "geom.h"
#include <vector>
#include<iostream>
#include<string>
#include <algorithm>
#include<cmath>
#include <set>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

template <typename T> int eps_sgn(T val) 
{
	if (eps_eq(val,T(0)))
	{
		return 0;
	}
	return (T(0) < val) - (val < T(0));
}


Vector vectprod(const Vector vector1, const Vector vector2)
{
	return { vector1.y() * vector2.z() - vector2.y() * vector1.z(),
			 vector1.z() * vector2.x() - vector2.z() * vector1.x(),
			 vector1.x() * vector2.y() - vector2.x() * vector1.y()};
}

double scalprod(const Vector vector1, const Vector vector2)
{
	return vector1.x() * vector2.x() + vector1.y() * vector2.y() + vector1.z() * vector2.z();
}
double sqrdist(const Point point1, const Point point2)
{
	return (pow(point2.x() - point1.x(), 2) + pow(point2.y() - point1.y(), 2) + pow(point2.z() - point1.z(), 2));
}
bool operator== (const Point point1, const Point point2)
{
	return eps_eq(point1.x(),point2.x()) && eps_eq(point1.y(),point2.y()) && eps_eq(point1.z(), point2.z());
}
bool operator== (const Segment segment1, const Segment segment2)
{
	return ((segment1.m_edges[0]==segment2.m_edges[0])&&(segment1.m_edges[1]==segment2.m_edges[1]))||
	       ((segment1.m_edges[0]==segment2.m_edges[1])&&(segment1.m_edges[1]==segment2.m_edges[0]));
}
Point operator*(const Point point, const double num)
{
	return { point.x() * num, point.y() * num, point.z() * num };
}
Point operator*(const double num, const Point point)
{
	return { point.x() * num, point.y() * num, point.z() * num };
}
Vector operator*(const double num, const Vector vector)
{
	return { vector.x() * num, vector.y() * num, vector.z() * num };
}
Vector operator*(const Vector vector, const double num)
{
	return { vector.x() * num, vector.y() * num, vector.z() * num };
}
bool operator==(const Vector vector1, const Vector vector2)
{
	return eps_eq(vector1.x(), vector2.x()) && eps_eq(vector1.y(), vector2.y()) && eps_eq(vector1.z(), vector2.z());
}
Vector operator+(const Vector vector1, const Vector vector2)
{
	return { vector1.x() + vector2.x(),vector1.y() + vector2.y() ,vector1.z() + vector2.z() };
}



	
int Triangle::actual() const
{
	if ((m_apexes[0] == m_apexes[1]) && (m_apexes[1] == m_apexes[2]))
	{
		return 0;
	}

	if (norm() == Vector{ 0,0,0 })
	{
		return 1;
	}

		return 2;

}
Segment Triangle::segm() const
{
	for (int i=0; i<3; i++)
	{
		if ((sqrdist(m_apexes[i], m_apexes[(i+1)%3])>=sqrdist(m_apexes[(i+1)%3], m_apexes[(i+2)%3]))&&(sqrdist(m_apexes[i], m_apexes[(i+1)%3])>=sqrdist(m_apexes[(i+2)%3], m_apexes[i])))
			{
				return Segment(m_apexes[i], m_apexes[(i+1)%3]);
			}
	}
}
Point Triangle::minCorner() const
{
	return {std::min(std::min(m_apexes[0].x(), m_apexes[1].x()), m_apexes[2].x()),
		    std::min(std::min(m_apexes[0].y(), m_apexes[1].y()), m_apexes[2].y()),
		    std::min(std::min(m_apexes[0].z(), m_apexes[1].z()), m_apexes[2].z())};
};
Point Triangle::maxCorner() const
{
	return {std::max(std::max(m_apexes[0].x(), m_apexes[1].x()), m_apexes[2].x()),
		    std::max(std::max(m_apexes[0].y(), m_apexes[1].y()), m_apexes[2].y()),
		    std::max(std::max(m_apexes[0].z(), m_apexes[1].z()), m_apexes[2].z())};
};
Point intsec(const Segment segment, const Triangle triangle)
{
	Vector temp = segment.m_edges[0]  + segment.dir() * scalprod({ segment.m_edges[0],triangle.m_apexes[0] }, triangle.norm())*(1/scalprod(segment.dir(), triangle.norm()));
	return { temp.x(), temp.y(), temp.z() };
}
std::istream& operator>> (std::istream& in, Triangle& triangle)
{
	double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	in >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
	triangle = Triangle{ x1, y1, z1, x2, y2, z2, x3, y3, z3 };
	return in;
}

bool compare(const Point point, const Segment segment)
{
	return (eps_eq(std::pow(scalprod(Vector{ segment.m_edges[0], point }, Vector{ segment.m_edges[1], point }), 2) , sqrdist(segment.m_edges[0], point) * sqrdist(segment.m_edges[1], point)) && (eps_sgn(scalprod(Vector{ segment.m_edges[0], point }, Vector{ segment.m_edges[1], point })) <= 0));
}
bool compare(const Segment segment1, const Segment segment2)
{
	if (eps_sgn(scalprod({ segment1.m_edges[0],segment2.m_edges[0] }, vectprod(segment1.dir(), segment2.dir()))) != 0)
	{
		return false;
	}
	Vector n1 = vectprod(segment1.dir(), vectprod(segment1.dir(), segment2.dir()));
	Vector n2 = vectprod(segment2.dir(), vectprod(segment1.dir(), segment2.dir()));
	int d11 = eps_sgn(scalprod(n2, { segment2.m_edges[0],segment1.m_edges[0] }));
	int d12 = eps_sgn(scalprod(n2, { segment2.m_edges[0],segment1.m_edges[1] }));
	int d21 = eps_sgn(scalprod(n1, { segment1.m_edges[0],segment2.m_edges[0] }));
	int d22 = eps_sgn(scalprod(n1, { segment1.m_edges[0],segment2.m_edges[1] }));
	return ((d11 != d12) && (d21 != d22)) ||
		((d11 == 0) && (d12 == 0) && ((compare(segment1.m_edges[0], segment2)) || (compare(segment1.m_edges[1], segment2)) || (compare(segment2.m_edges[0], segment1))));
}
bool compare(const Point point, const Triangle triangle)
{
	if (eps_sgn(scalprod({ triangle.m_apexes[0],point }, triangle.norm()))!=0)
	{
		return false;
	}
	return triangle.norm() ==
		Triangle{ point, triangle.m_apexes[0], triangle.m_apexes[1] }.norm() * eps_sgn(scalprod(triangle.norm(), Triangle{ point, triangle.m_apexes[0], triangle.m_apexes[1] }.norm())) +
		Triangle{ point, triangle.m_apexes[1], triangle.m_apexes[2] }.norm() * eps_sgn(scalprod(triangle.norm(), Triangle{ point, triangle.m_apexes[1], triangle.m_apexes[2] }.norm())) +
		Triangle{ point, triangle.m_apexes[2], triangle.m_apexes[0] }.norm() * eps_sgn(scalprod(triangle.norm(), Triangle{ point, triangle.m_apexes[2], triangle.m_apexes[0] }.norm()));
}
bool compare(const Segment segment, const Triangle triangle)
{
	int d1 = eps_sgn(scalprod({ triangle.m_apexes[0],segment.m_edges[0] }, triangle.norm()));
	int d2 = eps_sgn(scalprod({ triangle.m_apexes[0], segment.m_edges[1] }, triangle.norm()));
	if ((d1 + d2 == 2) || (d1 + d2 == -2))
	{
		return false;
	}
	else if ((d1 == 0) && (d2 == 0))
	{
		return compare(segment.m_edges[0], triangle) ||
			compare(segment.m_edges[1], triangle) ||
			compare(segment, { triangle.m_apexes[0],triangle.m_apexes[1] }) ||
			compare(segment, { triangle.m_apexes[1],triangle.m_apexes[2] }) ||
			compare(segment, { triangle.m_apexes[2],triangle.m_apexes[0] });
	}
	else if ((d1 == 0) && (d2 != 0))
	{
		if (compare(segment.m_edges[0], triangle))
		{
			return true;
		}
	}
	else if ((d1 != 0) && (d2 == 0))
	{
		if (compare(segment.m_edges[1], triangle))
		{
			return true;
		}
	}
	return compare(intsec(segment, triangle), triangle);
}
bool compare(const Triangle triangle1, const Triangle triangle2)
{
	if ((triangle1.actual() == 0) && (triangle2.actual() == 0))
	{
		return triangle1.m_apexes[0] == triangle2.m_apexes[0];
	}
	if ((triangle1.actual() == 0) && (triangle2.actual() == 1))
	{
		return (compare(triangle1.m_apexes[0], triangle2.segm()));
	}
	if ((triangle1.actual() == 1) && (triangle2.actual() == 0))
	{
		return (compare(triangle2.m_apexes[0], triangle1.segm()));
	}
	if ((triangle1.actual() == 1) && (triangle2.actual() == 1))
	{
		return (compare(triangle1.segm(), triangle2.segm()));
	}
	if ((triangle1.actual() == 0) && (triangle2.actual() == 2))
	{
		return (compare(triangle1.m_apexes[0], triangle2));
	}
	if ((triangle1.actual() == 2) && (triangle2.actual() == 0))
	{
		return (compare(triangle2.m_apexes[0], triangle1));
	}
	if ((triangle1.actual() == 1) && (triangle2.actual() == 2))
	{
		return (compare(triangle1.segm(), triangle2));
	}
	if ((triangle1.actual() == 2) && (triangle2.actual() == 1))
	{
		return (compare(triangle2.segm(), triangle1));
	}
	if ((triangle1.actual() == 2) && (triangle2.actual() == 2))
	{
		std::vector<int> d1(3);
		std::vector<int> d2(3);
		for (int i = 0; i < 3; i++)
		{
			d1[i] = eps_sgn(scalprod({ triangle2.m_apexes[0],triangle1.m_apexes[i] }, triangle2.norm()));
		}
		for (int i = 0; i < 3; i++)
		{
			d2[i] = eps_sgn(scalprod({ triangle1.m_apexes[0] , triangle2.m_apexes[i] }, triangle1.norm()));
		}
		if (((d1[0] == d1[1]) && (d1[1] == d1[2]) && (d1[0] != 0)) || ((d2[0] == d2[1]) && (d2[1] == d2[2]) && (d2[0] != 0)))
		{
			return false;
		}
		if ((d1[0] == 0) && (d1[1] == 0) && (d1[2] == 0))
		{
			return compare(triangle1.m_apexes[0], triangle2) || compare(triangle1.m_apexes[1], triangle2) || compare(triangle1.m_apexes[2], triangle2) ||
				compare(triangle2.m_apexes[0], triangle1) || compare(triangle2.m_apexes[1], triangle1) || compare(triangle2.m_apexes[2], triangle1);
		}
		//Having examined these 2 cases, we understand, the polygons intersect each other's planes, but do not lie on them
		Segment segment11, segment12, segment21, segment22;
		Point t11, t12, t21, t22;
		int k1{ 0 }, k2{ 0 };
		for (int i = 0; i < 3; i++)
		{
			//every apex of triangle 1 doesn't lie on the triangle 2's plane:
			if ((d1[i] == -d1[(i + 1) % 3]) && (d1[(i + 1) % 3] == d1[(i + 2) % 3]))
			{
				k1 = 2;
				segment11 = { triangle1.m_apexes[i], triangle1.m_apexes[(i + 1) % 3] };
				segment12 = { triangle1.m_apexes[(i + 2) % 3], triangle1.m_apexes[i] };
				t11 = intsec(segment11, triangle2);
				t12 = intsec(segment12, triangle2);
				break;
			}
			//only 2 apexes of triangle 1 lie on the plane of triangle 2:
			if ((d1[i] == 0) && (d1[(i + 1) % 3] == 0))
			{
				k1 = 2;
				t11 = triangle1.m_apexes[i];
				t12 = triangle1.m_apexes[(i + 1) % 3];
				break;
			}
			//only 1 apex of triangle 1 lie on the plane of triangle 2:
			if ((d1[i] == 0) && (d1[(i + 1) % 3] == -d1[(i + 2) % 3]))
			{
				k1 = 2;
				t11 = triangle1.m_apexes[i];
				segment11 = { triangle1.m_apexes[(i + 1) % 3] ,triangle1.m_apexes[(i + 2) % 3] };
				t12 = intsec(segment11, triangle2);
				break;
			}

			if ((d1[i] == 0) && (d1[(i + 1) % 3] == d1[(i + 2) % 3]))
			{
				k1 = 1;
				t11 = triangle1.m_apexes[i];
				break;
			}

		}
		for (int i = 0; i < 3; i++)
		{
			if ((d2[i] == -d2[(i + 1) % 3]) && (d2[(i + 1) % 3] == d2[(i + 2) % 3]))
			{
				k2 = 2;
				segment21 = { triangle2.m_apexes[i], triangle2.m_apexes[(i + 1) % 3] };
				segment22 = { triangle2.m_apexes[(i + 2) % 3], triangle2.m_apexes[i] };
				t21 = intsec(segment21, triangle1);
				t22 = intsec(segment22, triangle1);
				break;
			}

			if ((d2[i] == 0) && (d2[(i + 1) % 3] == 0))
			{
				k2 = 2;
				t21 = triangle2.m_apexes[i];
				t22 = triangle2.m_apexes[(i + 1) % 3];
				break;
			}

			if ((d2[i] == 0) && (d2[(i + 1) % 3] == -d2[(i + 2) % 3]))
			{
				k2 = 2;
				t21 = triangle2.m_apexes[i];
				segment21 = { triangle2.m_apexes[(i + 1) % 3] ,triangle2.m_apexes[(i + 2) % 3] };
				t22 = intsec(segment21, triangle1);
				break;
			}

			if ((d2[i] == 0) && (d2[(i + 1) % 3] == d2[(i + 2) % 3]))
			{
				k2 = 1;
				t21 = triangle2.m_apexes[i];
				break;
			}
		}
		if ((k1 == 1) && (k2 == 1))
		{
			return ((t21) == (t11));
		}
		else if ((k1 == 1) && (k2 == 2))
		{
			return compare(t11, Segment{ t21, t22});
			}
		else if ((k1 == 2) && (k2 == 1))
		{
			return compare(t21, Segment{ t11 , t12  });
			}
		else if ((k1 == 2) && (k2 == 2))
		{
			return compare(Segment{ t11, t12 }, Segment{ t21 , t22 });
		}
	}
}


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


void final()
{
	int N;
    std::cin >> N;
    std::vector<Triangle> triangles(N);
    std::vector<AABB> boundingBoxes(N);
    std::set<int> intersectingTriangles;
    
    using Box = bg::model::box<bg::model::point<double, 3, bg::cs::cartesian>>;
    using Value = std::pair<Box, int>;
    bgi::rtree<Value, bgi::quadratic<16>> rtree;
    
    for (int i = 0; i < N; ++i) {
        std::cin >> triangles[i];
        boundingBoxes[i] = {triangles[i].minCorner(), triangles[i].maxCorner()};
        Box box(
            bg::model::point<double, 3, bg::cs::cartesian>(boundingBoxes[i].minCorner.x(), boundingBoxes[i].minCorner.y(), boundingBoxes[i].minCorner.z()),
            bg::model::point<double, 3, bg::cs::cartesian>(boundingBoxes[i].maxCorner.x(), boundingBoxes[i].maxCorner.y(), boundingBoxes[i].maxCorner.z())
        );
        rtree.insert({box, i});
    }
    
    for (int i = 0; i < N; ++i) {
        Box queryBox(
            bg::model::point<double, 3, bg::cs::cartesian>(boundingBoxes[i].minCorner.x(), boundingBoxes[i].minCorner.y(), boundingBoxes[i].minCorner.z()),
            bg::model::point<double, 3, bg::cs::cartesian>(boundingBoxes[i].maxCorner.x(), boundingBoxes[i].maxCorner.y(), boundingBoxes[i].maxCorner.z())
        );
        
        std::vector<Value> candidates;
        rtree.query(bgi::intersects(queryBox), std::back_inserter(candidates));
        
        for (const auto& [box, j] : candidates) {
            if (i != j && compare(triangles[i], triangles[j])) {
                intersectingTriangles.insert(i);
                intersectingTriangles.insert(j);
            }
        }
    }
    
    for (int index : intersectingTriangles) {
        std::cout << index <<std::endl;
    }
}
