#ifndef LAST_STEP_INFO_INCLUDED
#define LAST_STEP_INFO_INCLUDED

#include "ray.hpp"
#include "edgestruct.hpp"
#include "propagation.hpp"
#include "typedef.hpp"

class LastStepInfo 
{
public:
    typedef ManifoldGraphT::Handle Handle;

public:
    // for DRAWING purposes
    int from_face, to_face;
    Interval* cur_iv;
    Handle cur_e;
    Vec2 s;
    Vec3 i0_pt, i1_pt;
    Ray<Vec2> ray0, ray1;
    Vec3 other_v;

    // modified intervals for 
    vector< pair<Interval*, Propagation> > modified_ivs;

public:
    LastStepInfo() 
    { from_face=to_face=-1; cur_iv=NULL; cur_e=NULL; }
    
    void reset() 
    { from_face = to_face = -1; cur_iv = NULL; cur_e = NULL; modified_ivs.clear(); }
};

#endif