#ifndef GEOM_H
#define GEOM_H
#include<iostream>
#include <vector>
#include <algorithm>
#include<cassert>
#include <iomanip>

template <typename T> bool eps_eq(T val1, T val2)
{
	T eps{0.000001};
	return std::abs(val2-val1)<eps;
}
template <typename T> int eps_sgn(T val);

class Point
{
private:
	std::vector<double> m_coords;
public:
	Point()=default;
	Point(const double x, const double y, const double z) : m_coords{ x,y,z } {}
	double x() const
	{
		return m_coords[0];
	}
	double y() const
	{
		return m_coords[1];
	}
	double z() const
	{
		return m_coords[2];
	}
};

class Vector
{
private:
	std::vector<double> m_coords;
public:
	Vector()=default;
	Vector(const double x, const double y, const double z) : m_coords{ x,y,z } {};
	Vector(const Point point1, const Point point2) : m_coords{ point2.x() - point1.x(),point2.y() - point1.y(),point2.z() - point1.z() } {};
	Vector(const Point point) : m_coords{ point.x(), point.y(), point.z() } {};
	double x() const
	{
		return m_coords[0];
	}
	double y() const
	{
		return m_coords[1];
	}
	double z() const
	{
		return m_coords[2];
	}
};

class Triangle;
Vector vectprod(const Vector vector1, const Vector vector2);

class Segment
{
private:
	std::vector<Point> m_edges;
public:
	Segment() = default;
	Segment(const Point edge1, const Point edge2) : m_edges{ std::vector<Point> {edge1, edge2} } {}
	Vector dir() const
	{
		return { m_edges[0], m_edges[1] };
	}
	friend bool operator== (const Segment segment1, const Segment segment2);
	friend bool compare(const Point point, const Segment segment);
	friend bool compare(const Segment segment1, const Segment segment2);
	friend bool compare(const Segment segment, const Triangle triangle);
	friend Point intsec(const Segment segment, const Triangle triangle);
};


class Triangle
{
private:
	std::vector <Point> m_apexes;

public:
	Triangle() {}
	Triangle(const Point apex1, const Point apex2, const Point apex3) : m_apexes{ std::vector <Point>{apex1, apex2, apex3} } {}
	Triangle(const double x1, const double y1, const double z1, const double x2, const double y2, const double z2, const double x3, const double y3, const double z3) :
		m_apexes{ std::vector<Point>{Point {x1, y1, z1},Point {x2, y2, z2},Point {x3, y3, z3}} } {}
	
	Vector norm() const
	{
		return  vectprod({ m_apexes[0], m_apexes[1] }, { m_apexes[1], m_apexes[2] });
	}
	int actual() const;
	Segment segm() const;
	Point minCorner() const;
	Point maxCorner() const;
	friend bool compare(const Triangle triangle1, const Triangle triangle2);
	friend bool compare(const Segment segment, const Triangle triangle);
	friend bool compare(const Point point, const Triangle triangle);
	friend Point intsec(const Segment segment, const Triangle triangle);
};


struct AABB {
    Point minCorner, maxCorner;
};


double scalprod(const Vector vector1, const Vector vector2);
double sqrdist(const Point point1, const Point point2);
bool operator== (const Point point1, const Point point2);
bool operator== (const Segment segment1, const Segment segment2);
Point operator*(const Point point, const double num);
Point operator*(const double num, const Point point);
Vector operator*(const double num, const Vector vector);
Vector operator*(const Vector vector, const double num);
bool operator==(const Vector vector1, const Vector vector2);
Vector operator+(const Vector vector1, const Vector vector2);

Point intsec(const Segment segment, const Triangle triangle);
std::istream& operator>> (std::istream& in, Triangle& triangle);


bool compare(const Point point, const Segment segment);
bool compare(const Segment segment1, const Segment segment2);
bool compare(const Point point, const Triangle triangle);
bool compare(const Segment segment, const Triangle triangle);
bool compare(const Triangle triangle1, const Triangle triangle2);
void final();
#endif
