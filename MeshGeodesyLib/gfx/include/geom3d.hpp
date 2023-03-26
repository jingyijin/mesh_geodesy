#ifndef GFXGEOM3D_INCLUDED
#define GFXGEOM3D_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  Handy 3D geometrical primitives

  $Id: geom3d.h 432 2004-11-02 22:55:41Z garland $

 ************************************************************************/

#include "vec3.hpp"

namespace gfx
{

//
// Computing properties of triangles
//

template<class Vec>
inline Vec triangle_raw_normal(const Vec& v1, const Vec& v2, const Vec& v3)
{
    return cross(v2-v1, v3-v1);
}

template<class Vec>
inline typename Vec::value_type
    triangle_area(const Vec& v1,const Vec& v2,const Vec& v3)
{
    return 0.5 * norm(triangle_raw_normal(v1, v2, v3));
}

template<class Vec>
inline Vec triangle_normal(const Vec& v1, const Vec& v2, const Vec& v3)
{
    Vec n = triangle_raw_normal(v1, v2, v3);
    unitize(n);
    return n;
}

template<class Vec, class Plane>
inline Plane triangle_plane(const Vec& v1, const Vec& v2, const Vec& v3)
{
    Vec n = triangle_normal(v1, v2, v3);
    return Plane(n, -(n*v1));
}

//
// Operations with axis-aligned bounding boxes
//

template<class Vec, class List>
void update_bbox(Vec& min, Vec& max, const List& items)
{
    typedef typename List::const_iterator iterator;

    for(iterator i=items.begin(); i!=items.end(); i++)
    {
    const Vec& v = *i;
    for(int j=0; j<Vec::dim(); j++)
    {
        if( v[j] < min[j] )  min[j] = v[j];	    
        if( v[j] > max[j] )  max[j] = v[j];
    }
    }
}

template<class Vec, class List>
void compute_bbox(Vec& min, Vec& max, const List& items)
{
    if( items.size()==0 )  min = max = 0;
    else                   min = max = items[0];

    update_bbox(min, max, items);
}

} // namespace gfx

// GFXGEOM3D_INCLUDED
#endif
