#ifndef TYPEDEF_INCLUDED
#define TYPEDEF_INCLUDED

/************************************************************************
 * File description: Typedefs
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "trimesh.hpp"
#include "manifold.hpp"
#include "interval.hpp"

typedef ManifoldGraph<TriMesh> ManifoldGraphT;
typedef ManifoldGraphT::Handle Handle;

typedef vector<Interval*> IntervalList;
typedef pair<double, double> Range;

#endif