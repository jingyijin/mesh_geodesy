#ifndef PATH_TRACER_INCLUDED
#define PATH_TRACER_INCLUDED

/************************************************************************
 * File description: Path tracer class to compute geodesic distances and paths
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "geomesh.hpp"

/**
 * @brief A class for tracing geodesic paths on a triangular mesh.
 */
class PathTracer
{
public:
    typedef vector<double> ScalarVector;
    typedef vector<ScalarVector> ScalarVectorVector;
    typedef GeoTriMesh::InterStruct InterStruct;
    typedef GeoTriMesh::HandleSet HandleSet;

    typedef pair<Handle, double> PathKnotPair;
    typedef vector<PathKnotPair> PathKnotVector;
    typedef vector<PathKnotVector> PathKnotVectorVector;

    enum State_type {Init_state, Tracing_ring_build_state, Tracing_ring_trace_state, 
        Tracing_path_state, Test_s_state, Test_sigma_state, Complete_state};

public:
    GeoTriMesh& m_mesh;       /**< The core mesh with attributes to calculate distance and paths. */
    PathKnotVectorVector& m_path_knot;    /**< The geodesic paths from a selected vertex. */

    State_type m_state;   /**< The state of the path tracer. */
    int m_next_vertex, m_target_vertex, m_vindex;   /**< The next vertex to be traced. */
    Interval* m_min_iv;   /**< The interval with the minimum distance. */
    Handle m_min_e;       /**< The edge with the minimum distance. */
    Handle m_cur_e;       /**< The current edge. */
    Handle m_prev_e;      /**< The previous edge. */
    InterStruct m_inter;  /**< The intersection structure. */
    Vec2 m_s;             /**< The pseudi start coordinate. */
    Vec2 m_start_coord;   /**< The start coordinate. */
    Vec2 m_inter_p2D;     /**< The intersection point in 2D. */
    int m_to_face;        /**< The face to be traced. */
    Ray<Vec2> m_ray;      /**< The ray to be traced. */
    map<Handle, bool> m_current_ring; /**< The current ring. */

public:
    /**
     * @brief Constructor for the PathTracer class.
     * @param m Reference to the input mesh.
     * @param path_knot Reference to the path knots.
     */
    PathTracer(GeoTriMesh& m, PathKnotVectorVector& path_knot);

    /**
     * @brief Clear the path data.
     */
    void clear_path();
    /**
     * @brief Trace a path.
     */
    void trace_path();
    /**
     * @brief Performs one step of path tracing using a finite state machine.
     * This method traces the path on the mesh by iterating through different states:
     * - Init_state: Initializes the tracing process.
     * - Tracing_ring_build_state: Builds the current ring.
     * - Tracing_ring_trace_state: Traces the ring to find the next vertex.
     * - Tracing_path_state: Traces the path between vertices.
     * - Test_s_state: Checks if the start and end coordinates are equal.
     * - Test_sigma_state: Checks if sigma is zero, indicating the end of the path.
     * - Complete_state: Marks the end of the current path tracing step.
     */
    void trace_once();

    void build_current_ring(int vid, Handle e, map<Handle, bool>& current_ring);
    void flag_current_ring(map<Handle, bool>& current_ring, Handle e);
    Interval* trace_ray_once(int to_face, Handle cur_e, Ray<Vec2>& ray, InterStruct& inter);
    Interval* get_one_interval(map<Handle, bool>& current_ring, Vec2& start_coord);
    Interval* get_interval_at(Handle& eh, double inter);
    int get_interval_end_vertex(Interval* iv, Vec2& start_coord);
    int get_next_vertex(Interval* iv);

    void store_path_knot(int source_v, int knot_index, Handle eh=NULL);
    void store_path_knot(int source_v, Handle e, double interval);
};


#endif