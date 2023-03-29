//////////////////////////////////////////////////////////////////////////////////////
//	
//////////////////////////////////////////////////////////////////////////////////////

#ifndef GMATH_INCLUDED
#define GMATH_INCLUDED

#include "ray.hpp"
#include "vec2.hpp"
#include "vec3.hpp"
#include "mat4.hpp"

#include <set>
#include <vector>

#define BETWEEN(a,b,c) (a>b && a<c)

using namespace std;

class GMath
{
public:
    struct Line : pair<Vec3, Vec3> {
        Line(Vec3& p0, Vec3& p1) { first = p0; second = p1; }
        Vec3 p() { return first; }
        Vec3 d() { Vec3 dir = second-first; unitize(dir); return dir; }
    };
    struct Segment : pair<Vec3,Vec3> {
        Line supporting_line() { return Line(first, second); }
    };

public:
    static Vec3 rotate_point(const Vec3& moving, double angle, const Vec3& axis, const Vec3& offset);
    static Vec3 rotate_direction(const Vec3& moving, double angle, const Vec3& axis);
    static double compute_signed_angle_deg(const Vec3 v0, const Vec3& v1, const Vec3& ref);
    static double compute_dihedral_angle(const Vec3 v0, const Vec3& v1);

    static bool within_range(double interv, const pair<double, double>& range);
    static bool collinear(const Vec3& p, const Vec3& q, const Vec3& r);
    static bool are_ordered_along_line(const Vec2& p, const Vec2& q, const Vec2& r);
    static bool are_ordered_along_line(const Vec3& p, const Vec3& q, const Vec3& r);

    static bool vector_equal(Vec3& v0, Vec3& v1);
    static bool is_null(Vec3& v);
    static double dot(Vec3& p, Vec3& q, Vec3& r);
    static bool is_acute_angle(Vec3& u, Vec3& v);
    static bool is_acute_angle(Vec3& p, Vec3& q, Vec3& r);
    static short sign_of_determinant2x2(const double& a00, const double& a01, const double& a10, const double& a11);
    static short compare(const double& n1, const double& n2);

    static double area(const Vec2& p0, const Vec2& p1, const Vec2& p2);
    static Vec3 barycenter(const Vec2& p0, const Vec2& p1, const Vec2& p2, const Vec2& pc);
    static bool all_positive(const Vec3& v);
    static Vec2 convert_to_local(const Vec3& ref, const Vec3& global);
    static Vec3 convert_to_global(const Vec3& ref, const Vec2& local);

    // Solve quadratic function
    static int solve(const double a, const double b, const double c, vector<double>& res);
    
    // return the conversion multiplicative from degree to radians
    static double degree_to_radians() 
    {	
        return 0.017453292;		
    }

    // return the conversion multiplicative from radians to degree
    static double radians_to_degree() 
    {
        return 57.2957795131;	
    }

    static double angle_between(const Vec3& v0, const Vec3& v1)
    {
        return asin(norm(v0 ^ v1));
    }

    static double unit_dotp(const Vec3& v0, const Vec3& v1)
    {
        Vec3 temp_v0 = v0, temp_v1 = v1;
        unitize(temp_v0);
        unitize(temp_v1);

        return (temp_v0 * temp_v1);
    }

    static double distance_to_plane(const Vec3& n, double d, const Vec3& point)
    {
        return ((n*point)+d);
    }

    static double sign(double value)
    {
        return (value < 0) ? -1 : 1;
    }
};

#endif