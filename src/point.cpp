#include "point.hpp"

#include <cassert>
#include <cmath>
#include <cfloat>




Point::Point()
{
  for (int i = 0; i < N_COORD; ++i)
    _coord[i] = 0.;
}



Point::Point(const double coordinates[])
{
  for (int i = 0; i < N_COORD; ++i)
    _coord[i] = coordinates[i];
}



Point::Point(const std::vector<double> &coordinates)
{
  assert(N_COORD == (int)coordinates.size());

  for (int i = 0; i < N_COORD; ++i)
    _coord[i] = coordinates[i];
}



Point::Point(double x_coord,
             double y_coord,
             double z_coord)
{
  _coord[0] = x_coord;
  _coord[1] = y_coord;
  _coord[2] = z_coord;
}



Point::Point(const Point &p)
{
  for (int i = 0; i < N_COORD; ++i)
    _coord[i] = p._coord[i];
}



Point& Point::operator =(const Point &p)
{
  for (int i = 0; i < N_COORD; ++i)
    _coord[i] = p._coord[i];
  return *this;
}



double Point::operator ()(int number) const
{
  assert(number < N_COORD);
  return coord(number);
}



double Point::coord(int number) const
{
  assert(number < N_COORD);
  return _coord[number];
}



double Point::x() const
{
  return _coord[0];
}



double Point::y() const
{
  assert(N_COORD > 1);
  return _coord[1];
}



double Point::z() const
{
  assert(N_COORD > 2);
  return _coord[2];
}



double& Point::operator ()(int number)
{
  assert(number < N_COORD);
  return coord(number);
}



double& Point::coord(int number)
{
  assert(number < N_COORD);
  return _coord[number];
}



double& Point::x()
{
  return _coord[0];
}



double& Point::y()
{
  assert(N_COORD > 1);
  return _coord[1];
}



double& Point::z()
{
  assert(N_COORD > 2);
  return _coord[2];
}



void Point::coord(int number, double value)
{
  assert(number < N_COORD);
  _coord[number] = value;
}



Point& Point::operator +=(const Point &p)
{
  for (int i = 0; i < N_COORD; ++i)
    _coord[i] += p._coord[i];
  return *this;
}



Point& Point::operator *=(double d)
{
  for (int i = 0; i < N_COORD; ++i)
    _coord[i] *= d;
  return *this;
}



Point& Point::operator /=(double d)
{
  assert(fabs(d) > DBL_EPSILON*DBL_EPSILON);

  for (int i = 0; i < N_COORD; ++i)
    _coord[i] /= d;

  return *this;
}



Point operator +(const Point &p1, const Point &p2)
{
  Point res;
  for (int i = 0; i < Point::N_COORD; ++i)
    res._coord[i] = p1._coord[i] + p2._coord[i];
  return res;
}



Point operator -(const Point &p1, const Point &p2)
{
  Point res;
  for (int i = 0; i < Point::N_COORD; ++i)
    res._coord[i] = p1._coord[i] - p2._coord[i];
  return res;
}



std::ostream& operator <<(std::ostream &os, const Point &p)
{
  os << "(";
  for (int i = 0; i < Point::N_COORD; ++i)
    os << p._coord[i] << (i == Point::N_COORD - 1 ? "" : ",");
  os << ")";
  return os;
}



bool Point::compare_by_x(const Point &a, const Point &b)
{
  return (a._coord[0] < b._coord[0]);
}



bool Point::compare_by_y(const Point &a, const Point &b)
{
  return (a._coord[1] < b._coord[1]);
}



bool Point::compare_by_z(const Point &a, const Point &b)
{
  return (a._coord[2] < b._coord[2]);
}



bool Point::operator ()(const Point &a, const Point &b) const
{
  return compare_by_norm(a, b);
}



bool Point::operator <(const Point &p) const
{
  return compare_by_norm(*this, p);
}



bool Point::compare_by_norm(const Point &a, const Point &b)
{
  const Point c = a - b;
  double norm_a = 0., norm_b = 0., norm_c = 0.;
  for (int i = 0; i < N_COORD; ++i)
  {
    norm_a += a._coord[i] * a._coord[i];
    norm_b += b._coord[i] * b._coord[i];
    norm_c += c._coord[i] * c._coord[i];
  }
  norm_a = sqrt(norm_a);
  norm_b = sqrt(norm_b);
  norm_c = sqrt(norm_c);

  const double tol = DBL_EPSILON;
  return (norm_c > tol && norm_a < norm_b);
}


