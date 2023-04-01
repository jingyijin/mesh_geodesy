#ifndef CELL_TYPES_H_INCLUDED 
#define CELL_TYPES_H_INCLUDED

/************************************************************************
 * File description: Cell types for mesh generation
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include <set>
#include <vector>
#include <cstddef>

typedef unsigned int iid_t;
typedef std::pair<iid_t, iid_t> idpair_t;
typedef std::set<iid_t> idset_t;
typedef std::vector<iid_t> idlist_t;

struct Polygon : public idlist_t
{
    Polygon(size_t N) : idlist_t(N) {}
    Polygon(const idlist_t& p) : idlist_t(p) {}

    int dim() const { return 2; }
};

typedef idset_t VertSet;
typedef idset_t EdgeSet;

#endif