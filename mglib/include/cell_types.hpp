#ifndef CELL_TYPES_H_INCLUDED 
#define CELL_TYPES_H_INCLUDED

#include <set>
#include <vector>

typedef std::pair<id_t, id_t> idpair_t;
typedef std::set<id_t> idset_t;
typedef std::vector<id_t> idlist_t;

struct Polygon : public idlist_t
{
	Polygon(size_t N) : idlist_t(N) {}
	Polygon(const idlist_t& p) : idlist_t(p) {}

	int dim() const { return 2; }
};

typedef idset_t VertSet;
typedef idset_t EdgeSet;

#endif