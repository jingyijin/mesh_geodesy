#ifndef RAY_INCLUDED 
#define RAY_INCLUDED

#include "vec3.hpp"

template <class Vec>
class Ray
{
public:
    enum { Infinite_intersection, One_intersection, No_intersection } intersection_type;

public:
    Vec p;
    Vec d;
    double t;	// the distance from the origin to the point used to specify the ray's direction

public:
    Ray() {}

    Ray(const Vec& source, const Vec& target)
    {
        makeRay(source, target);
    }

    void makeRay(const Vec& source, const Vec& target)
    {
        p = source;
        d = target - source;
        t = norm(d);
        unitize(d);
    }

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

    bool has_infinite_intersection()const { return intersection_type == Infinite_intersection; }
    bool has_one_intersection()		const { return intersection_type == One_intersection; }
    bool has_no_intersection()		const { return intersection_type == No_intersection; }

    double sign(double value)
    {
        return (value < 0) ? -1 : 1;
    }
};

#endif