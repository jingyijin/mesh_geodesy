#ifndef TYPEDEF_INCLUDED
#define TYPEDEF_INCLUDED

#include "trimesh.hpp"
#include "manifold.hpp"
#include "interval.hpp"

typedef ManifoldGraph<TriMesh> ManifoldGraphT;
typedef ManifoldGraphT::Handle Handle;

typedef vector<Interval*> IntervalList;
typedef pair<double, double> Range;

#endif