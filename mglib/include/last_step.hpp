#ifndef LAST_STEP_INFO_INCLUDED
#define LAST_STEP_INFO_INCLUDED

/************************************************************************
 * File description: Last step information for mesh generation
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/


#include "ray.hpp"
#include "propagation.hpp"


/**
 * @brief A class to hold information about the last propagation step.
 * This class is used for DRAWING purposes and to store modified intervals during propagation.
*/
class LastStepInfo 
{
public:
    typedef ManifoldGraphT::Handle Handle;

public:
    int from_face, to_face; /**< Indexes of the faces before and after the last propagation step */
    Interval *cur_iv;       /**< The current interval */
    Handle cur_e;           /**< Handle to the current halfedge */
    Vec2 s;                 /**< Psudo source point */
    Vec3 i0_pt, i1_pt;      /**< Intersection points */
    Ray<Vec2> ray0, ray1;   /**< Rays */
    Vec3 other_v;           /**< Other vector */

    // modified intervals for 
    vector< pair<Interval*, Propagation> > modified_ivs;    /**< Modified intervals */

public:
    /**
     * @brief Default constructor.
     */
    LastStepInfo() 
        : cur_iv(nullptr), cur_e(NULL), from_face(-1), to_face(-1) 
    {}

    /**
     * @brief Default destructor.
     */
    ~LastStepInfo() {
        if (cur_iv) delete cur_iv;
    }
    /**
     * @brief Resets the object to its default state.
     */
    void reset() 
    { from_face = to_face = -1; cur_e = NULL; modified_ivs.clear(); }
};

#endif