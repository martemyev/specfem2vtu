#ifndef FEMPLUS_POINT_HPP
#define FEMPLUS_POINT_HPP

#include <iostream>
#include <vector>




/**
 * @brief A 3D point
 *
 * A class to describe a point in 3-dimensional space. It's used for
 * 2-dimensional points as well, and in this case one of the coordinates is 0
 * (usually it's z-coordinate).
 *
 * @author M. Artemyev
 * @date 2014, 2015
 */
class Point
{
public:

  /**
   * @brief Number of coordinates describing a point
   *
   * The number of Cartesian coordinates, that describe the point. Here we
   * always use 3 coordinates to describe a point.
   */
  static const int N_COORD = 3;

  /**
   * @brief Default constructor creates a point at origin
   *
   * Default constructor. Coordinates are initialized by 0.
   */
  Point();

  /**
   * @brief Costructor based on array of coordinates
   *
   * Constructor with parameter. Coordinates are initialized by a given array of
   * numbers.
   * @param coordinates - array of point coordinates
   */
  Point(const double coordinates[]);

  /**
   * @brief Constructor based on a vector of coordinates
   *
   * Constructor with parameter. Coordinates are initialized by a std vector of
   * values.
   * @param coordinates - a vector of point coordinates
   */
  Point(const std::vector<double> &coordinates);

  /**
   * @brief Constructor based on concrete values of coordinates
   *
   * Constructor with parameters. Coordinates are initialized by numbers.
   * @param x_coord - x-coordinate of the point
   * @param y_coord - y-coordinate of the point
   * @param z_coord - z-coordinate of the point
   */
  Point(double x_coord,
        double y_coord = 0,
        double z_coord = 0);

  /**
   * @brief Copy constructor
   *
   * Copy constructor
   *  @param p - a point from which the copy is made
   */
  Point(const Point &p);

  /**
   * @brief Copy assignment operator
   *
   * Copy assignment operator
   * @param p - a point from which the copy is made
   * @return newly created (copied) point
   */
  Point& operator =(const Point &p);

  /**
   * @brief point_a += point_b
   *
   * Add a point to this one: point_a = point_a + point_b.
   * The points are added by coordinates.
   * @param p - a point to be added
   */
  Point& operator +=(const Point &p);

  /**
   * @brief point *= scalar
   *
   * Multiply the point by some number.
   * In this case all coordinates are multiplied by this number.
   * It is a scaling.
   * @param d - scaling factor
   */
  Point& operator *=(double d);

  /**
   * @brief point /= scalar
   *
   * Divide the point by some number.
   * In this case all coordinates are divided by this number.
   * It is a scaling.
   * @param d - denominator
   */
  Point& operator /=(double d);

  /**
   * @brief |point1| < |point2|
   *
   * Comparison between two points based on their norms (i.e. distance from the
   * origin). This function simply calls compare_by_norm(const Point &a, const Point &b).
   * @param a - one point for comparison
   * @param b - another point for comparison
   * @return true if 'b' point is farther from the origin than 'a' point
   */
  bool operator ()(const Point &a, const Point &b) const;

  /**
   * @brief |this_point| < |another_point|
   *
   * Comparison between this and another points based on their norms (i.e.
   * distance from the origin). This function simply calls
   * compare_by_norm(const Point &a, const Point &b).
   * @param p - another point for comparison
   * @return true if 'p' point is farther from the origin than 'this' point
   */
  bool operator <(const Point &p) const;

  /**
   * @brief Get a coordinate
   *
   * Get a coordinate of the point. Simply calls coord(int number) const.
   * @param number - a serial number of the coordinate of interest
   * @return a value of the coordinate
   */
  double operator ()(int number) const;

  /**
   * @brief Get/set a coordinate
   *
   * Get/set a coordinate of the point. Simply calls coord(int number).
   * @param number - a serial number of the coordinate of interest
   * @return a reference to the coordinate so it can be changed
   */
  double& operator ()(int number);

  /**
   * @brief point3 = point1 + point2
   *
   * Add one point to another one: p = p1 + p2.
   * The points are added by coordinates.
   * @param p1 - one point (left operand)
   * @param p2 - another point (right operand)
   * @return p = p1 + p2
   */
  friend Point operator +(const Point &p1, const Point &p2);

  /**
   * @brief point3 = point1 - point2
   *
   * Subtract one point from another one: p = p1 - p2.
   * The points are substracted by coordinates.
   * @param p1 - one point (left operand)
   * @param p2 - another point (right operand)
   * @return p = p1 - p2
   */
  friend Point operator -(const Point &p1, const Point &p2);

  /**
   * @brief Output to a stream
   *
   * Output by '<<' to some stream
   * @param os - output stream (e.g. std::cout)
   * @param p - the point itself
   * @return changed output stream
   */
  friend std::ostream& operator <<(std::ostream &os, const Point &p);

  /**
   * @brief Get a coordinate
   *
   * Get a coordinate of the point
   * @param number - the serial number of coordinate [0, N_COORD)
   * @return a value of the coordinate (its copy)
   */
  double coord(int number) const;

  /**
   * @brief Get an x-coordinate
   *
   * Get an x-coordinate of the point
   * @return x-coordinate of the point
   */
  double x() const;

  /**
   * @brief Get a y-coordinate
   *
   * Get a y-coordinate of the point
   * @return y-coordinate of the point
   */
  double y() const;

  /**
   * @brief Get a z-coordinate
   *
   * Get a z-coordinate of the point
   * @return z-coordinate of the point
   */
  double z() const;

  /**
   * @brief Get/set a coordinate
   *
   * Get/set a coordinate of the point
   * @param number - the serial number of coordinate [0, N_COORD)
   * @return a reference to a corresponding coordinate (the coordinate can be
   *         changed then)
   */
  double& coord(int number);

  /**
   * @brief Get/set an x-coordinate
   *
   * Get/set an x-coordinate of the point
   * @return a reference to an x-coordinate of the point
   */
  double& x();

  /**
   * @brief Get/set a y-coordinate
   *
   * Get/set a y-coordinate of the point
   * @return a reference to a y-coordinate of the point
   */
  double& y();

  /**
   * @brief Get/set a z-coordinate
   *
   * Get/set a z-coordinate of the point
   * @return a reference to a z-coordinate of the point
   */
  double& z();

  /**
   * @brief Set a coordinate
   *
   * Set the value of specific coordinate
   * @param number - the number of coordinate that we want to set
   * @param value - new value of coordinate
   */
  void coord(int number, double value);

  /**
   * @brief point1.x < point2.x
   *
   * The comparison between points which is based on their x-coordinates.
   * The bigger x-coordinates, the bigger the point.
   * This function is used, for example, to sort the vector of points
   * according to their x-coordinates.
   * @param a - one point
   * @param b - another point
   * @return true if b-point has bigger x-coordinate
   */
  static bool compare_by_x(const Point &a, const Point &b);

  /**
   * @brief point1.y < point2.y
   *
   * The comparison between points which is based on their y-coordinates.
   * The bigger y-coordinates, the bigger the point.
   * This function is used, for example, to sort the vector of points
   * according to their y-coordinates.
   * @param a - one point
   * @param b - another point
   * @return true if b-point has bigger y-coordinate
   */
  static bool compare_by_y(const Point &a, const Point &b);

  /**
   * @brief point1.z < point2.z
   *
   * The comparison between points which is based on their z-coordinates.
   * The bigger z-coordinates, the bigger the point.
   * This function is used, for example, to sort the vector of points
   * according to their z-coordinates.
   * @param a - one point
   * @param b - another point
   * @return true if b-point has bigger z-coordinate
   */
  static bool compare_by_z(const Point &a, const Point &b);

  /**
   * @brief |point1| < |point2|
   *
   * The comparison between points which is based on the norm of the
   * vector between the point and the origin (0, 0, 0).
   * The farther from origin, the bigger the point.
   * This function is used, for example, to sort the vector of points
   * according to their distance from the origin.
   * @param a - one point
   * @param b - another point
   * @return true if b-point has bigger norm
   */
  static bool compare_by_norm(const Point &a, const Point &b);


private: // ================= PRIVATE =======================

  /**
   * Cartesian coordinates of the point
   */
  double _coord[N_COORD];
};



#endif // FEMPLUS_POINT_HPP
