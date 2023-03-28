#ifndef PROPAGATION_H_INCLUDED
#define PROPAGATION_H_INCLUDED

#include "edgestruct.hpp"

class Propagation
{
public:
    EdgeStruct deleted;
    EdgeStruct altered;

public:
    Propagation() {}

    Propagation(EdgeStruct& al) 
    { altered.resize(al.size()); copy(al.begin(), al.end(), altered.begin()); }

    void clear() 
    { deleted.clear(); altered.clear(); }

    void clear_deleted()
    { deleted.clear(); }

    void clear_altered()
    { altered.clear(); }

    void insert_altered(IntervalList& ivs) 
    { altered.insert_intervals(ivs); }

    void insert_altered(Interval* iv) 
    { altered.insert_interval(iv); }

    void insert_deleted(Interval* iv)
    { deleted.insert_interval(iv); }
};

#endif