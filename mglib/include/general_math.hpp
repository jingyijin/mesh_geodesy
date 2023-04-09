#ifndef GMATH_INCLUDED
#define GMATH_INCLUDED

/************************************************************************
 * File description: Math functions for geometry
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "vec2.hpp"
#include "vec3.hpp"
#include "mat4.hpp"

#include <set>
#include <vector>
#include <glog/logging.h>

#define BETWEEN(a,b,c) (a>b && a<c)

const double FEQ_EPS_LOOSE = 1e-3;

inline bool FLT(double a, double b, double e=FEQ_EPS)  { return (b-a)>e;}
inline bool FGT(double a, double b, double e=FEQ_EPS)  { return (a-b)>e;}
inline bool FLT_LOOSE(double a, double b, double e=FEQ_EPS_LOOSE)  { return (b-a)>e;}
inline bool FGT_LOOSE(double a, double b, double e=FEQ_EPS_LOOSE)  { return (a-b)>e;}

extern double current_time, prev_time;

using namespace std;

/**
* @class GMath
* @brief A collection of static methods for geometric computations.
*/
class GMath
{
public:
    /**
    * @brief A line in 3D space, defined by two points.
    */
    struct Line : pair<Vec3, Vec3> {
        Line(Vec3& p0, Vec3& p1) { first = p0; second = p1; }
        Vec3 p() { return first; }
        Vec3 d() { Vec3 dir = second-first; unitize(dir); return dir; }
    };

    /**
    * @brief A line segment in 3D space, defined by two points.
    */
    struct Segment : pair<Vec3,Vec3> {
        Line supporting_line() { return Line(first, second); }
    };

public:
    /**
    * @brief Rotate a point in 3D space around a given axis by a given angle.
    * @param moving The point to be rotated.
    * @param angle The rotation angle in radians.
    * @param axis The axis of rotation.
    * @param offset The offset from the origin of the rotation axis.
    * @return The rotated point.
    */
    static Vec3 rotate_point(const Vec3& moving, double angle, const Vec3& axis, const Vec3& offset);
    /**
    * @brief Rotate a direction vector in 3D space around a given axis by a given angle.
    * @param moving The direction vector to be rotated.
    * @param angle The rotation angle in radians.
    * @param axis The axis of rotation.
    * @return The rotated direction vector.
    */    
    static Vec3 rotate_direction(const Vec3& moving, double angle, const Vec3& axis);
    /**
    * @brief Compute the signed angle in degrees between two vectors, with respect to a reference vector.
    * @param v0 The first vector.
    * @param v1 The second vector.
    * @param ref The reference vector.
    * @return The signed angle between the two vectors in degrees.
    */
    static double compute_signed_angle_deg(const Vec3 v0, const Vec3& v1, const Vec3& ref);
    /**
    * @brief Compute the dihedral angle between two vectors.
    * @param v0 The first vector.
    * @param v1 The second vector.
    * @return The dihedral angle between the two vectors.
    */
    static double compute_dihedral_angle(const Vec3 v0, const Vec3& v1);

    /**
    * @brief Check if a value is within a given range.
    * @param interv The value to check.
    * @param range The range to check against.
    * @return True if the value is within the range, false otherwise.
    */
    static bool within_range(double interv, const pair<double, double>& range);
    /**
    * @brief Check if three points are collinear.
    * @param p The first point.
    * @param q The second point.
    * @param r The third point.
    * @return True if the points are
    */
    static bool collinear(const Vec3& p, const Vec3& q, const Vec3& r);
    /**
    * Determines if points p, q, and r are ordered along a line.
    * @param p The first point
    * @param q The second point
    * @param r The third point
    * @return true if the points are ordered along a line, false otherwise
    */    
    static bool are_ordered_along_line(const Vec2& p, const Vec2& q, const Vec2& r);
    /**
    * Determines if points p, q, and r are ordered along a line.
    * @param p The first point
    * @param q The second point
    * @param r The third point
    * @return true if the points are ordered along a line, false otherwise
    */    
    static bool are_ordered_along_line(const Vec3& p, const Vec3& q, const Vec3& r);

    /**
    * Determines if two 3D vectors are equal.
    * @param v0 The first 3D vector
    * @param v1 The second 3D vector
    * @return true if the vectors are equal, false otherwise
    */
    static inline bool vector_equal(Vec3& v0, Vec3& v1);
    /**
    * @brief Check if a 3D vector is a null vector.
    * @param v A reference to the 3D vector to be checked.
    * @return True if the 3D vector is a null vector, false otherwise.
    */
    static inline bool is_null(Vec3& v);
    /**
    * @brief Calculate the dot product between two 3D vectors, given a reference point.
    * @param p A reference to the point.
    * @param q A reference to the first 3D vector.
    * @param r A reference to the second 3D vector.
    * @return The dot product between the two 3D vectors.
    */
    static inline double dot(Vec3& p, Vec3& q, Vec3& r);
    /**
    * @brief Check if two 3D vectors form an acute angle.
    * @param u A reference to the first 3D vector.
    * @param v A reference to the second 3D vector.
    * @return True if the two 3D vectors form an acute angle, false otherwise.
    */    
    static inline bool is_acute_angle(Vec3& u, Vec3& v);
    /**
    * @brief Check if three points form an acute angle.
    * @param p A reference to the first point.
    * @param q A reference to the second point.
    * @param r A reference to the third point.
    * @return True if the three points form an acute angle, false otherwise.
    */    
    static inline bool is_acute_angle(Vec3& p, Vec3& q, Vec3& r);
    /**
    * Determines the sign of a 2x2 matrix determinant.
    * @param a00 The 0,0 entry of the matrix
    * @param a01 The 0,1 entry of the matrix
    * @param a10 The 1,0 entry of the matrix
    * @param a11 The 1,1 entry of the matrix
    * @return -1 if the determinant is negative, 0 if the determinant is zero, 1 if the determinant is positive
    */        
    static short sign_of_determinant2x2(const double& a00, const double& a01, const double& a10, const double& a11);
    /**
    * Compares two double values.
    * @param n1 The first double value
    * @param n2 The second double value
    * @return -1 if n1 < n2, 0 if n1 = n2, 1 if n1 > n2
    */
    static short compare(const double& n1, const double& n2);

    /**
    * @brief Computes the area of a triangle given its vertices in 2D space.
    * @param p0 first vertex of the triangle
    * @param p1 second vertex of the triangle
    * @param p2 third vertex of the triangle
    * @return double the area of the triangle
    */
    static double area(const Vec2& p0, const Vec2& p1, const Vec2& p2);
    /**
    * @brief Computes the barycenter of a triangle and a point.
    * @param p0 first vertex of the triangle
    * @param p1 second vertex of the triangle
    * @param p2 third vertex of the triangle
    * @param pc the point to compute the barycenter with
    * @return Vec3 the barycenter of the triangle and point
    */    
    static Vec3 barycenter(const Vec2& p0, const Vec2& p1, const Vec2& p2, const Vec2& pc);
    /**
    * Checks if all the components of a Vec3 vector are positive or zero.
    * @param v The vector to check.
    * @return True if all the components of the vector are positive or zero, false otherwise.
    */
    static bool all_positive(const Vec3& v);
    /**
    * Converts a global point to a local point with respect to a reference vector.
    * @param ref The reference vector.
    * @param global The global point to convert.
    * @return A Vec2 representing the local point.
    */
    static Vec2 convert_to_local(const Vec3& ref, const Vec3& global);
    /**
    * Converts a local point to a global point with respect to a reference vector.
    * @param ref The reference vector.
    * @param local The local point to convert.
    * @return A Vec3 representing the global point.
    */
    static Vec3 convert_to_global(const Vec3& ref, const Vec2& local);

    /**
    * @brief Solves a quadratic function of the form ax^2 + bx + c = 0.
    * @param a coefficient of x^2
    * @param b coefficient of x
    * @param c constant term
    * @param res vector to store the roots
    * @return int number of roots
    */
    static int solve(const double a, const double b, const double c, vector<double>& res);
    
    /**
    * @brief Returns the conversion multiplicative factor from degrees to radians.
    * @return The factor for converting degrees to radians (approximately 0.017453292).
    */
   static inline double degree_to_radians() 
    {
        return 0.017453292;
    }

    /**
    * @brief Returns the conversion multiplicative factor from radians to degrees.
    * @return The factor for converting radians to degrees (approximately 57.2957795131).
    */
    static inline double radians_to_degree() 
    {
        return 57.2957795131;	
    }

    /**
    * @brief Calculates the angle in radians between two vectors.
    * @param v0 The first vector.
    * @param v1 The second vector.
    * @return The angle in radians between the two vectors.
    */
    static inline double angle_between(const Vec3& v0, const Vec3& v1)
    {
        return asin(norm(v0 ^ v1));
    }

    /**
    * @brief Calculates the dot product between two unit vectors.
    * @param v0 The first unit vector.
    * @param v1 The second unit vector.
    * @return The dot product between the two unit vectors.
    */
    static inline double unit_dotp(const Vec3& v0, const Vec3& v1)
    {
        Vec3 temp_v0 = v0, temp_v1 = v1;
        unitize(temp_v0);
        unitize(temp_v1);

        return (temp_v0 * temp_v1);
    }

    /**
    * @brief Calculates the signed distance from a point to a plane.
    * @param n The normal vector of the plane.
    * @param d The distance from the plane to the origin.
    * @param point The point to calculate the distance to.
    * @return The signed distance from the point to the plane.
    */
    static inline double distance_to_plane(const Vec3& n, double d, const Vec3& point)
    {
        return ((n*point)+d);
    }

    /**
    * @brief Returns the sign of a given value.
    * @param value The value to determine the sign of.
    * @return -1 if the value is negative, 1 if the value is positive.
    */
    static inline double sign(double value)
    {
        return (value < 0) ? -1 : 1;
    }
};

template <class T>
inline void print_vector(const vector<T>& vec)
{
    cout << "vector size: " << vec.size() << endl;
    for (typename vector<T>::const_iterator iter = vec.begin(); iter != vec.end(); iter++)
        cout << *iter << endl;
}

template <class T>
inline void copy_vector(vector<T>& from, vector<T>& to)
{
    to.clear();
    to.resize(from.size());
    copy(from.begin(), from.end(), to.begin());
}

template <class T>
inline void attach_vector(vector<T>& first, vector<T>& second)
{
    for (typename vector<T>::iterator lit = second.begin(); lit != second.end(); lit++)
        first.push_back(*lit);
}

inline void reset_time()
{
    prev_time = get_cpu_time();
}

inline double get_used_time()
{
    current_time = get_cpu_time();
    double diff = current_time - prev_time;
    prev_time = current_time;
    return diff;
}

#endif