//////////////////////////////////////////////////////////////////////////////////////
//	Auxiliary file that implements methods of class GMath (General Math)			//
//////////////////////////////////////////////////////////////////////////////////////

#include "general_math.hpp"

Vec3 GMath::rotate_point(const Vec3& moving, double angle, const Vec3& axis, const Vec3& offset)
{
    Mat4 M = rotation_matrix_deg(angle, axis);

    Vec4 v(moving-offset, 0.f);
    Vec4 vnew = M*v;

    return Vec3(vnew[0]+offset[0], vnew[1]+offset[1], vnew[2]+offset[2]);
}

/* rotate a 3D vector by angle about axis */
Vec3 GMath::rotate_direction(const Vec3& moving, double angle, const Vec3& axis)
{
    Mat4 M = rotation_matrix_deg(angle, axis);
    Vec4 v = Vec4f(moving, 0.0);
    Vec4 vnew = M*v;

    return Vec3(vnew[0], vnew[1], vnew[2]);
}

double GMath::compute_signed_angle_deg(const Vec3 v0, const Vec3& v1, const Vec3& ref)
{
    Vec3 cross_n = v0 ^ v1;
    double angle = 0, v_norm2 = v0*v1;
    if (v_norm2 < 1)
        angle = GMath::radians_to_degree() * acos(v0*v1);
    if ((cross_n*ref) < 0.0) angle = -angle;
    return angle;
}

double GMath::compute_dihedral_angle(const Vec3 v0, const Vec3& v1)
{
    Vec3 cross_n = v0 ^ v1;
    return GMath::radians_to_degree() * asin(norm(cross_n));
}

bool GMath::within_range(double interv, const pair<double, double>& range)
{
    return !(FLT(interv, range.first) || FGT(interv, range.second));
}
/*
bool GMath::collinear(const Vec3& p, const Vec3& q, const Vec3& r)
{
    // following CGAL implementation
    Vec3 dp = p-r;
    Vec3 dq = q-r;
    if (GMath::sign_of_determinant2x2(dp[0], dq[0], dp[1], dq[1]) == 0)
        return false;

    return GMath::sign_of_determinant2x2(dp[0], dq[0], dp[2], dq[2]) == 0 && 
           GMath::sign_of_determinant2x2(dp[1], dq[1], dp[2], dq[2]) == 0;
}
*/
/*
bool GMath::collinear(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    Vec3 d1 = p2-p1;
    Vec3 d2 = p3-p1;
    Vec3 c = d1^d2;

    return FEQ(c*c, 0.0);
}
*/
bool GMath::collinear(const Vec3& a, const Vec3& b, const Vec3& c)
{
    double det = a * (b ^ c);

    if (FEQ(det, 0))
    {
        // test further
        // when one of three points is linearly dependent of the others
        // the determinant is 0, but this does not imply that they are collinear
        Vec3 ba = b - a;
        Vec3 ca = c - a;
        unitize(ba);
        unitize(ca);
        return FEQ(fabs(ba * ca), 1.0);
    }

    return false;
}

bool GMath::are_ordered_along_line(const Vec2& p, const Vec2& q, const Vec2& r)
{
    if (FLT(p[0],q[0])) return !FLT(r[0],q[0]);
    if (FLT(q[0],p[0])) return !FLT(q[0],r[0]);
    if (FLT(p[1],q[1])) return !FLT(r[1],q[1]);
    if (FLT(q[1],p[1])) return !FLT(q[1],r[1]);
    return true; // p==q
}

bool GMath::are_ordered_along_line(const Vec3& p, const Vec3& q, const Vec3& r)
{
    if (FLT(p[0],q[0])) return !FLT(r[0],q[0]);
    if (FLT(q[0],p[0])) return !FLT(q[0],r[0]);
    if (FLT(p[1],q[1])) return !FLT(r[1],q[1]);
    if (FLT(q[1],p[1])) return !FLT(q[1],r[1]);
    if (FLT(p[2],q[2])) return !FLT(r[2],q[2]);
    if (FLT(q[2],p[2])) return !FLT(q[2],r[2]);
    return true; // p==q
}

short GMath::sign_of_determinant2x2(const double& a00, const double& a01, const double& a10, const double& a11)
{
    return compare( a00*a11, a10*a01);
}

short GMath::compare(const double& n1, const double& n2)
{
    return (n1 < n2) ? -1 : (n2 < n1) ? 1 : 0;
}

bool GMath::vector_equal(Vec3& v0, Vec3& v1)
{
    return FEQ(norm2(v0-v1), 0.0);
}

bool GMath::is_null(Vec3& v)
{
    return FEQ(norm2(v), 0.0);
}

double GMath::dot(Vec3& p, Vec3& q, Vec3& r)
{
  return  ( (p[0] - q[0]) * (r[0] - q[0])
        + (p[1] - q[1]) * (r[1] - q[1])
        + (p[2] - q[2]) * (r[2] - q[2]) );
}

bool GMath::is_acute_angle(Vec3& u, Vec3& v)
{
    return (u*v) > double(0) ;
}

bool GMath::is_acute_angle(Vec3& p, Vec3& q, Vec3& r)
{
    return GMath::dot(p, q, r) > double(0) ;
}

int GMath::solve(const double a, const double b, const double c, vector<double>& res) 
{
    int num_roots = 0;
    double delta = b*b-4*a*c;
    double half = -b/(2*a);
    if (delta < 0) {
        num_roots = 0;
    } else if (FEQ(delta, 0.0)) {
        num_roots = 1;
        res.push_back(half);
    } else {
        num_roots = 2;
        delta = sqrt(delta);
        res.push_back(half + (delta/(2*a)));
        res.push_back(half - (delta/(2*a)));
    }
    return num_roots;
}

Vec3 GMath::barycenter(const Vec2& p0, const Vec2& p1, const Vec2& p2, const Vec2& pc) 
{
    Vec3 bary;
    double area012 = GMath::area(p0, p1, p2);
    bary[0] = GMath::area(pc, p1, p2)/area012;
    bary[1] = GMath::area(pc, p2, p0)/area012;
    bary[2] = 1.0 - bary[0] - bary[1];

    return bary;
}

double GMath::area(const Vec2& p0, const Vec2& p1, const Vec2& p2)
{
    return ((p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1]))/2.0;
} 

bool GMath::all_positive(const Vec3& v)
{
    for (int i=0; i<3; i++)
        if (FLT(v[i], 0.0)) 
            return false;
    return true;
}

Vec2 GMath::convert_to_local(const Vec3& ref, const Vec3& global)
{
    Vec3 u, v(0,0,1);
    if (FEQ(norm2(ref-v), 0.0))
        v = Vec3(0,1,0);
    u = ref ^ v;
    v = ref ^ u;
    unitize(u);
    unitize(v);

    return Vec2(global*u, global* v);
}

Vec3 GMath::convert_to_global(const Vec3& ref, const Vec2& local)
{
    Vec3 u, v(0,0,1);
    if (FEQ(norm2(ref-v), 0.0))
        v = Vec3(0,1,0);
    u = ref ^ v;
    v = ref ^ u;
    unitize(u);
    unitize(v);

    return local[0]*u + local[1]*v;
}
