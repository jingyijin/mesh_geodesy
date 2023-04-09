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

/**
 * @brief A class to compute geodesic distances and paths on a triangular mesh.
 */
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
    GeoTriMesh *mesh;               /**< The core mesh with fundamental structures to compute geodesic. */

    // geodesic distances and paths
    ScalarVector distances;         /**< The geodesic distances from a selected vertex. */
    KnotVectorVector paths;         /**< The geodesic paths from a selected vertex. */

    ScalarVector scalar_field;      /**< The scalar field on the mesh. */

public:
    /**
     * @brief Constructor for the MeshGeodesy class.
     * @param m A pointer to the input mesh.
     */
    MeshGeodesy(GeoTriMesh *m);
    /**
     * @brief Destructor for the MeshGeodesy class.
     */
    ~MeshGeodesy();

    /**
     * @brief Clear all computed data.
     */
    void clear();
    /**
     * @brief Clear only computed distances.
     */
    void clear_distances();

    /**
     * @brief Compute geodesic distances from a selected vertex.
     * @param selected_v The index of the selected vertex.
     */
    void compute_distances(int selected_v);
    /**
     * @brief Sort mesh faces by their geodesic distance.
     */
    void sort_faces_by_distance();

    /**
     * @brief Save computed geodesic distances to a file.
     * @param filename The output file name.
     */
    void save_distances(const string& filename);
    /**
     * @brief Load geodesic distances from a file.
     * @param filename The input file name.
     */
    void load_distances(const string& filename);

    /**
     * @brief Print a knot vector for debugging purposes.
     * @param kn The input knot vector.
     */
    void print_knot_vector(const KnotVector& kn);
};

#endif