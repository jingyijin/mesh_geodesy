#ifndef MESH_GEODESY_INCLUDED
#define MESH_GEODESY_INCLUDED

/************************************************************************
 * File description: MeshGeodesy class to compute geodesic distances and paths
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "geomesh.hpp"

class MeshGeodesy
{
public:
    typedef vector< pair<int, double> > ValueConstraintMap;
    typedef vector< pair<int, Vec3> > GradientConstraintMap;
    typedef vector<double> ScalarVector;
    typedef vector<Vec3> VectorVector;

    typedef pair<Handle, double> KnotPair;
    typedef vector<KnotPair> KnotVector;
    typedef vector<KnotVector> KnotVectorVector;

public:
    GeoTriMesh *mesh;

    // geodesic distances and paths
    ScalarVector distances;
    KnotVectorVector paths;

    ScalarVector scalar_field;

public:
    MeshGeodesy(GeoTriMesh *m);
    ~MeshGeodesy();

    void clear();
    void clear_distances();

    void compute_distances(int selected_v);
    void sort_faces_by_distance();

    void save_distances(const string& filename);
    void load_distances(const string& filename);

    // in/output methods
    void print_knot_vector(const KnotVector& kn);
};

#endif