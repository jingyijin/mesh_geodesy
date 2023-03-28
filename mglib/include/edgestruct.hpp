#ifndef EDGESTRUCT_H_INCLUDED
#define EDGESTRUCT_H_INCLUDED

#include <list>
#include <assert.h>
#include <float.h>

#include "typedef.hpp"

class EdgeStruct : public vector<Interval*>
{
public:
    typedef Interval::Handle Handle;

    struct cmp_range
    {
        bool operator()(const Range& a, const Range& b) const
        {
            if (a.first == b.first)
                return a.second < b.second;
            return a.first < b.first;
        }
    };
    typedef set<Range, cmp_range> RangeSet;
    
public:
    bool modified;

public:
    EdgeStruct()
    { modified = false;	}
    EdgeStruct(Interval* iv)
    { push_back(iv); modified = false; }

    IntervalList::const_iterator get_const_begin() const { return begin(); }
    IntervalList::const_iterator get_const_end() const { return end(); }

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

    void insert_intervals(IntervalList& ivs) 
    {
        for (IntervalList::iterator iit = ivs.begin(); iit != ivs.end(); iit++)
            insert_interval(*iit);
    }

    void replace_interval_list(IntervalList& ivs)
    {
        clear();
        insert_intervals(ivs);
    }

    double get_distance_end0() const
    {
        if (empty())	return DBL_MAX;
        return front()->get_d0_dist();
    }

    double get_distance_end1() const
    {
        if (empty())	return DBL_MAX;
        return back()->get_d1_dist();
    }

    double get_smallest_b0() { return front()->get_b0(); }
    double get_largest_b1() { return back()->get_b1(); }

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

    // assertion methods
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

    bool has_flipped_interval()
    {
        bool flipped = false;
        for (IntervalList::iterator iit = begin(); iit != end() && !flipped; iit++)
            flipped = FGT((*iit)->get_b0(), (*iit)->get_b1());
        return flipped;
    }

    bool has_out_range_interval(double el)
    {
        if (empty()) return false;
        bool out_range = false;
        out_range = FGT_LOOSE(0.0, front()->get_b0());
        if (size() > 1 && !out_range)
            out_range = FGT_LOOSE(back()->get_b1(), el);
        return out_range;
    }

    bool has_flipped_edge() 
    {
        bool flipped = false;
        for (IntervalList::iterator iit = begin(); iit != end() && !flipped; iit++)
            flipped = (*iit)->handle()->Dest() < (*iit)->handle()->Org();
        return flipped;
    }

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


inline std::ostream& operator<<(std::ostream& out, const EdgeStruct& e)
{
    int i=0;
    for (IntervalList::const_iterator iit = e.get_const_begin(); iit != e.get_const_end(); iit++, i++)
    {
        out << "iv " << i << ": " << *(*iit) << endl;
    }
    return out;
}

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