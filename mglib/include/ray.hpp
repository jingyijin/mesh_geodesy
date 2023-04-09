#ifndef RAY_INCLUDED 
#define RAY_INCLUDED

/************************************************************************
 * File description: Ray class for flattened calculation of the wavefront
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "vec3.hpp"
#include "general_math.hpp"

/**
 * @brief A template class for a ray in n-dimensional space.
 * @tparam Vec A vector class representing points in n-dimensional space.
 */
template <class Vec>
class Ray
{
public:
    enum { Infinite_intersection, One_intersection, No_intersection } intersection_type;    /**< The type of intersection. */

public:
    Vec p;      /**< The origin of the ray.*/
    Vec d;      /**< The direction of the ray. */
    double t;   /**< The length of the ray. */

public:
    /**
     * @brief Default constructor.
     */
    Ray() {}
    /**
     * @brief Constructs a ray from source to target.
     * @param source The source point.
     * @param target The target point.
     */
    Ray(const Vec& source, const Vec& target)
    {
        makeRay(source, target);
    }

    /**
     * @brief Constructs a ray from source to target.
     * @param source The source point.
     * @param target The target point.
     */
    void makeRay(const Vec& source, const Vec& target)
    {
        p = source;
        d = target - source;
        t = norm(d);
        unitize(d);
    }

    /**
     * @brief Computes the intersection of the ray with a plane.
     * @param normal The normal vector of the plane.
     * @param D The distance from the origin to the plane.
     * @param[out] point The intersection point.
     * @return The distance from the ray's origin to the intersection point.
     */
    double planeIntersect(const Vec& normal, const double D, Vec& point)
    {
        double num = normal*p + D;
        double den = normal*d;
        if (FEQ(den, 0.0)) {
            if (FEQ(num, 0.0))
                intersection_type = Infinite_intersection;
            else 
                intersection_type = No_intersection;
            return 0;
        }

        intersection_type = One_intersection;
        double t = -(num)/(den);
        if (t < 0.0)
            return 0.0;

        t *= -sign(den);
        point = p + t*d;
        return t;
    }

    /**
     * @brief Computes the intersection of the ray with a line segment.
     * @param A The starting point of the line segment.
     * @param B The ending point of the line segment.
     * @param[out] point The intersection point.
     * @return The distance from the ray's origin to the intersection point.
     */
    double segment_intersect(const Vec& A, const Vec& B, Vec& point)
    {
        Vec C = p;
        Vec D = p+d;
        double r_nom = (A[1]-C[1])*(D[0]-C[0])-(A[0]-C[0])*(D[1]-C[1]);
        double r_den = (B[0]-A[0])*(D[1]-C[1])-(B[1]-A[1])*(D[0]-C[0]);

        if (FEQ(r_den, 0.0) && FEQ(r_nom, 0.0)) {
            intersection_type = Infinite_intersection;
            return 0.0;
        }

        double s_nom = (A[1]-C[1])*(B[0]-A[0])-(A[0]-C[0])*(B[1]-A[1]);
        double s_den = (B[0]-A[0])*(D[1]-C[1])-(B[1]-A[1])*(D[0]-C[0]);
        double r = r_nom/r_den;
        double s = s_nom/s_den;
        double c = (B[1]-A[1])*B[0] - (B[0]-A[0])*B[1];
        double side = (B[1]-A[1])*p[0] - (B[0]-A[0])*p[1] - c;
        if (!FLT(r, 0.0) && !FGT(r, 1.0) && !FLT(s, 0.0) && !FGT(side, 0)) {
            intersection_type = One_intersection;
            point = p + s*d;
            return s;
        }

        Vec BA = B-A; 
        unitize(BA);
        if (FEQ(fabs(BA*d), 1.0) && FEQ(side, 0.0)) {
            intersection_type = Infinite_intersection;
            return 0.0;
        }

        intersection_type = No_intersection;
        return 0.0;
    }

    /**
     * @brief Checks if the ray has an infinite intersection.
     * @return true if the ray has an infinite intersection, false otherwise.
     */
    inline bool has_infinite_intersection()const { return intersection_type == Infinite_intersection; }
    /**
     * @brief Checks if the ray has one intersection.
     * @return true if the ray has one intersection, false otherwise.
     */
    inline bool has_one_intersection()     const { return intersection_type == One_intersection; }
    /**
     * @brief Checks if the ray has no intersection.
     * @return true if the ray has no intersection, false otherwise.
     */
    inline bool has_no_intersection()      const { return intersection_type == No_intersection; }

    /**
     * @brief Returns the sign of the value.
     * @param value The input value.
     * @return -1 if the value is negative, 1 otherwise.
     */
    inline double sign(double value)
    {
        return (value < 0) ? -1 : 1;
    }
};

#endif