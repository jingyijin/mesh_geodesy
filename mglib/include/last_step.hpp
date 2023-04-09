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
    int m_from_face, m_to_face; /**< Indexes of the faces before and after the last propagation step */
    Interval *m_cur_iv;         /**< The current interval */
    Handle m_cur_e;             /**< Handle to the current halfedge */
    Vec2 m_s;                   /**< Psudo source point */
    Vec3 m_i0_pt, m_i1_pt;      /**< Intersection points */
    Ray<Vec2> m_ray0, m_ray1;   /**< Rays */
    Vec3 m_other_v;             /**< Other vector */

    // modified intervals for 
    vector< pair<Interval*, Propagation>> m_modified_ivs;    /**< Modified intervals */

public:
    /**
     * @brief Default constructor.
     */
    LastStepInfo() 
        : m_cur_iv(nullptr), m_cur_e(NULL), m_from_face(-1), m_to_face(-1) 
    {}

    /**
     * @brief Default destructor.
     */
    ~LastStepInfo() {
        if (m_cur_iv) delete m_cur_iv;
    }
    /**
     * @brief Resets the object to its default state.
     */
    void reset() 
    { m_from_face = m_to_face = -1; m_cur_e = NULL; m_modified_ivs.clear(); }
};

#endif