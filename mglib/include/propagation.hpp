#ifndef PROPAGATION_H_INCLUDED
#define PROPAGATION_H_INCLUDED

/************************************************************************
 * File description: Propagation class for the wavefront algorithm
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "edgestruct.hpp"

/**
 * @brief Propagation class for the wavefront algorithm.
 */
class Propagation
{
public:
    EdgeStruct m_deleted; /**< The edges that have been deleted during the propagation. */
    EdgeStruct m_altered; /**< The edges that have been altered during the propagation. */

public:
    /**
     * @brief Default constructor.
     */
    Propagation() {}
    /**
     * @brief Constructs a propagation object with an initial set of altered edges.
     * @param al The initial set of altered edges.
     */
    Propagation(EdgeStruct& al) 
    { m_altered.resize(al.size()); copy(al.begin(), al.end(), m_altered.begin()); }

    /**
     * @brief Clears both the deleted and altered edge sets.
     */
    inline void clear() 
    { m_deleted.clear(); m_altered.clear(); }
    /**
     * @brief Clears the deleted edge set.
     */
    inline void clear_deleted()
    { m_deleted.clear(); }
    /**
     * @brief Clears the altered edge set.
     */
    inline void clear_altered()
    { m_altered.clear(); }

    /**
     * @brief Inserts a list of intervals into the altered edge set.
     * @param ivs The list of intervals to insert.
     */
    void insert_altered(IntervalList& ivs) 
    { m_altered.insert_intervals(ivs); }
    /**
     * @brief Inserts a single interval into the altered edge set.
     * @param iv The interval to insert.
     */
    void insert_altered(Interval* iv) 
    { m_altered.insert_interval(iv); }
    /**
     * @brief Inserts a single interval into the deleted edge set.
     * @param iv The interval to insert.
     */
    void insert_deleted(Interval* iv)
    { m_deleted.insert_interval(iv); }
};

#endif