#ifndef EDGESTRUCT_H_INCLUDED
#define EDGESTRUCT_H_INCLUDED

/************************************************************************
 * File description: Edge structure for mesh generation
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include <list>
#include <assert.h>
#include <float.h>

#include "trimesh.hpp"
#include "manifold.hpp"
#include "interval.hpp"

typedef ManifoldGraph<TriMesh> ManifoldGraphT;
typedef ManifoldGraphT::Handle Handle;

typedef vector<Interval*> IntervalList;
typedef pair<double, double> Range;


/**
* @brief Class representing an edge with its intervals.
* This class represents an edge of a mesh with its intervals (represented as Interval objects).
* It provides methods for manipulating the intervals and checking for errors in the structure.
*/
class EdgeStruct : public vector<Interval*>
{
public:
    typedef Interval::Handle Handle;
    /**
    * @brief Functor used for comparing Range objects by their starting point and ending point.
    */
    struct cmp_range
    {
        bool operator()(const Range& a, const Range& b) const
        {
            if (a.first == b.first)
                return a.second < b.second;
            return a.first < b.first;
        }
    };
    /**
    * @brief Set of Range objects used for checking for uncovered parts of the edge.
    */
    typedef set<Range, cmp_range> RangeSet;
    
public:
    bool modified; /**< Flag indicating if the edge has been modified since it was last checked. */

public:
    /**
    * @brief Default constructor for EdgeStruct.
    */
    EdgeStruct()
    { modified = false;	}
    /**
     * @brief Constructor for EdgeStruct with one initial Interval.
     * 
     * @param iv The Interval to add to the EdgeStruct.
     */
    EdgeStruct(Interval* iv)
    { push_back(iv); modified = false; }
    /**
     * @brief Destructor for EdgeStruct.
     */    
    ~EdgeStruct() { 
        clear();
    }
    /**
     * @brief Returns a constant iterator to the beginning of the IntervalList.
     * @return A constant iterator to the beginning of the IntervalList.
     */
    IntervalList::const_iterator get_const_begin() const { return begin(); }
    /**
     * @brief Returns a constant iterator to the end of the IntervalList.
     * @return A constant iterator to the end of the IntervalList.
     */
    IntervalList::const_iterator get_const_end() const { return end(); }

    /**
     * @brief Inserts an Interval into the EdgeStruct.
     * Inserts an Interval into the EdgeStruct in order by the b0 coordinate of the Interval.
     * @param iv The Interval to insert into the EdgeStruct.
     */
    void insert_interval(Interval* iv) 
    {
        if (empty()) push_back(iv);
        else {
            bool inserted = false;
            IntervalList::iterator iit = this->begin();
            for ( ; iit != this->end() && !inserted; iit++) 
                if ((*iit)->get_b0() > iv->get_b0())
                {
                    insert(iit, iv);
                    inserted = true;
                }
            if (!inserted) push_back(iv);
        }
    }

    /**
     * @brief Inserts a list of Intervals into the EdgeStruct.
     * Inserts a list of Intervals into the EdgeStruct by calling insert_interval() on each Interval.
     * @param ivs The list of Intervals to insert into the EdgeStruct.
     */
    void insert_intervals(IntervalList& ivs) 
    {
        for (IntervalList::iterator iit = ivs.begin(); iit != ivs.end(); iit++)
            insert_interval(*iit);
    }

    /**
     * @brief Replaces the list of Intervals in the EdgeStruct.
     * Replaces the list of Intervals in the EdgeStruct with a new list of Intervals.
     * @param ivs The new list of Intervals to replace the old list with.
     */
    void replace_interval_list(IntervalList& ivs)
    {
        clear();
        insert_intervals(ivs);
    }

    /**
    * @brief Get the distance from the start of the edge to the closest point on the mesh for the start vertex
    * @return double The distance from the start of the edge to the closest point on the mesh for the start vertex
    */
    double get_distance_end0() const
    {
        if (empty())	return DBL_MAX;
        return front()->get_d0_dist();
    }

    /**
    * @brief Get the distance from the end of the edge to the closest point on the mesh for the end vertex
    * @return double The distance from the end of the edge to the closest point on the mesh for the end vertex
    */
    double get_distance_end1() const
    {
        if (empty())	return DBL_MAX;
        return back()->get_d1_dist();
    }

    /**
    * @brief Get the smallest b0 value of all intervals in the EdgeStruct
    * @return double The smallest b0 value of all intervals in the EdgeStruct
    */
    double get_smallest_b0() { return front()->get_b0(); }
    /**
    * @brief Get the largest b1 value of all intervals in the EdgeStruct
    * @return double The largest b1 value of all intervals in the EdgeStruct
    */    
    double get_largest_b1() { return back()->get_b1(); }

    /**
    * @brief Get the interval containing a specified point
    * @param inter The point of interest
    * @return Interval* The interval containing the specified point or the last interval in the EdgeStruct if the point is outside all intervals
    */
    Interval* get_interval(double inter)
    {
        for (IntervalList::iterator iit = begin(); iit != end(); iit++)
        {
            Interval* iv = (*iit);
            if (!FLT(iv->get_b1(), inter) && !FGT(iv->get_b0(), inter))
                return iv;
        }
        return back();
    }

    /**
    * @brief Check if there is any interval in the EdgeStruct with a null handle
    * @return true If there is any interval in the EdgeStruct with a null handle
    * @return false If there is no interval in the EdgeStruct with a null handle
    */
    bool has_null_element() 
    {
        bool has_null_element = false;
        for (IntervalList::iterator iit = begin(); iit != end() && !has_null_element; iit++)
        {
            if ((*iit)->handle() == NULL)
                has_null_element = true;
        }
        return has_null_element;
    }

    /**
    * @brief Check if the intervals in the EdgeStruct are out of order
    * @return true If the intervals in the EdgeStruct are out of order
    * @return false If the intervals in the EdgeStruct are in order
    */
    bool is_out_of_order()
    {
        if (size() < 2)
            return false;

        bool out_of_order = false;
        Interval* iv = front();
        double prev_b0 = iv->get_b0();
        for (IntervalList::const_iterator iit = begin()+1; iit != end() && !out_of_order; iit++)
        {
            iv = *iit;
            double current_b0 = iv->get_b0();
            out_of_order = FGT(prev_b0, current_b0);
            prev_b0 = current_b0;
        }
        return out_of_order;
    }

    /**
    * @brief Check if there is any overlapping interval in the EdgeStruct
    * @return true If there is any overlapping interval in the EdgeStruct
    * @return false If there is no overlapping interval in the EdgeStruct
    */
    bool has_overlapping_interval()
    {
        if (empty())
            return false;
        
        bool overlapping = false;
        Interval* prev = front();
        for (IntervalList::iterator iit = begin()+1; iit != end() && !overlapping; iit++)
            overlapping = (*iit)->does_intersect(prev);
        return overlapping;
    }


    /**
    * Check if the edge structure has any flipped interval.
    * @return true if the edge structure has a flipped interval, false otherwise.
    */
    bool has_flipped_interval()
    {
        bool flipped = false;
        for (IntervalList::iterator iit = begin(); iit != end() && !flipped; iit++)
            flipped = FGT((*iit)->get_b0(), (*iit)->get_b1());
        return flipped;
    }

    /**
    * Check if the edge structure has any interval outside a given range.
    * @param el The upper limit of the range.
    * @return true if the edge structure has an interval outside the given range, false otherwise.
    */
    bool has_out_range_interval(double el)
    {
        if (empty()) return false;
        bool out_range = false;
        out_range = FGT_LOOSE(0.0, front()->get_b0());
        if (size() > 1 && !out_range)
            out_range = FGT_LOOSE(back()->get_b1(), el);
        return out_range;
    }

    /**
    * Check if the edge structure has any flipped edge.
    * @return true if the edge structure has a flipped edge, false otherwise.
    */
    bool has_flipped_edge() 
    {
        bool flipped = false;
        for (IntervalList::iterator iit = begin(); iit != end() && !flipped; iit++)
            flipped = (*iit)->handle()->Dest() < (*iit)->handle()->Org();
        return flipped;
    }

    /**
    * Check if the edge structure covers a given range and store the uncovered ranges.
    * @param covering A reference to a RangeSet object where uncovered ranges are stored.
    * @param e_length The length of the edge structure.
    */
    void check_covering(RangeSet& covering, double e_length)
    {
        if (!empty())
        {
            double init = 0;
            for (IntervalList::iterator iit = begin(); iit != end(); iit++)
            {
                Interval* iv = *iit;
                if (FGT(iv->get_b0(), init))
                    covering.insert(make_pair(init, iv->get_b0()));
                init = iv->get_b1();
            }

            // do for the covering at the end
            Interval* iv = back();
            if (FLT(iv->get_b1(), e_length))
                covering.insert(make_pair(iv->get_b1(), e_length));
        }
    }
};


/**
* @brief Overloaded operator to print EdgeStruct object
* @param out output stream
* @param e EdgeStruct object to be printed
* @return output stream
*/
inline std::ostream& operator<<(std::ostream& out, const EdgeStruct& e)
{
    int i=0;
    for (IntervalList::const_iterator iit = e.get_const_begin(); iit != e.get_const_end(); iit++, i++)
    {
        out << "iv " << i << ": " << *(*iit) << endl;
    }
    return out;
}

/**
* @brief Overloaded operator to print a vector of Interval pointers
* @param out output stream
* @param ivs vector of Interval pointers to be printed
* @return output stream
*/
inline std::ostream& operator<<(std::ostream& out, const vector<Interval*>& ivs)
{ 
    int i=0;
    out << "vector size: " << ivs.size() << endl;
    for (vector<Interval*>::const_iterator iit = ivs.begin(); iit != ivs.end(); iit++, i++)
    {
        out << "iv " << i << ": " << *(*iit) << endl;
    }
    return out; 
}

#endif