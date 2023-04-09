#ifndef GEO_SMF_INCLUDED
#define GEO_SMF_INCLUDED

/************************************************************************
 * File description: GeoTriMesh class for geodesic computation. It holds
 * the key mesh data structure and the geodesic computation.
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "interval.hpp"
#include "manifold.hpp"
#include "heap.hpp"
#include "last_step.hpp"

#include <map>

typedef ManifoldGraph<TriMesh> ManifoldGraphT;
typedef Interval::Handle Handle;

/**
 * @brief A class for computing geodesic distances and paths on a triangular mesh.
 */
class GeoTriMesh : public TriMesh
{
public:
    typedef set<Handle> HandleSet;
    typedef EdgeStruct::RangeSet RangeSet;
    typedef vector<double> ScalarVector;
    typedef pair<int, int> IntPair;

    /**
     * @brief A comparator for comparing edge handles.
     */
    struct cmp_handle
    {
        bool operator()(const Handle a, const Handle b) const
        {
            if (a->Org() == b->Org())
                return a->Dest() < b->Dest();
            return a->Org() < b->Org();
        }
    };
    /**
     * @brief A comparator for comparing integer pairs.
     */
    struct cmp_intpair
    {
        bool operator()(const IntPair &a, const IntPair &b) const
        {
            if (a.first == b.first)
                return a.second < b.second;
            return a.first < b.first;
        }
    };
    /**
     * @brief A struct for storing information about an intersection point between two edges.
     */
    struct InterStruct
    {
        Handle e;     /**< The edge containing the intersection point. */
        double inter; /**< The distance along the edge where the intersection occurs. */
        double t;     /**< The normalized distance of the intersection point along the edge. */

        InterStruct()
        {
            e = NULL;
            inter = 0.;
            t = 0.;
        }
        /**
         * @brief Check whether the intersection point coincides with the end of an interval.
         *
         * @param e_value The endpoint value.
         * @return true if the intersection point is the same as the endpoint, false otherwise.
         */
        bool is_interp_end(double e_value) const { return FEQ(inter, e_value); }
    };

    typedef map<Handle, EdgeStruct, cmp_handle> EdgeStructMap;   /**< A map of edge handles to EdgeStructs. */
    typedef map<IntPair, RangeSet, cmp_intpair> EdgeCoveringMap; /**< A map of integer pairs to RangeSets. */

    typedef pair<Handle, Interval *> IntervalPair; /**< A pair of an edge handle and an Interval. */

    typedef pair<Handle, double> PathKnotPair;       /**< A pair of an edge handle and a distance along the edge. */
    typedef vector<PathKnotPair> PathKnotVector;     /**< A vector of PathKnotPairs. */
    typedef vector<PathKnotVector> KnotVectorVector; /**< A vector of PathKnotVectors. */

public:
    ManifoldGraphT *graph; /**< The manifold graph of the mesh. */

    // ******************************
    // for geodesic computation
    MxHeap i_queue;            /**< The heap for storing intervals. */
    EdgeStructMap edge_map;    /**< The map of edge handles to EdgeStructs. */
    ScalarVector geo_distance; /**< The vector of geodesic distances. */
    KnotVectorVector geo_path; /**< The vector of geodesic paths. */

    int prop_step_size;        /**< The propagation step size. */
    int num_prop_edge;         /**< The number of edges to propagate. */
    int stop_prop_number;      /**< The number of edges to stop propagation. */
    bool is_first_propagation; /**< Whether it is the first propagation. */
    double ERROR_TOLERANCE;    /**< The error tolerance for the propagation. */

    LastStepInfo last_step; /**< The last step information. */

    // for checking purpose
    EdgeCoveringMap edge_covering;

public:
    /**
    * @brief Constructor for GeoTriMesh.
    * @param m A pointer to the TriMesh object to be used as a base.
    This constructor initializes the GeoTriMesh by calling the parent constructor and then calling the initialize() method.
    */
    GeoTriMesh(TriMesh *m);
    /**
     * @brief Destructor for GeoTriMesh.
     * This destructor resets the distance properties and deletes the ManifoldGraphT object if it exists.
     */
    ~GeoTriMesh();

    /**
     * @brief Initializes the GeoTriMesh object.
     * This method creates a ManifoldGraphT object and initializes the properties related to geodesic computations,
     * such as the propagation step size, the number of propagated edges, the geodesic distances vector, and more.
     */
    void initialize();
    /**
     * @brief Resets the distance properties of the GeoTriMesh object.
     * This method resets the properties related to geodesic distances and paths, such as the edge map, the MxHeap,
     * the geodesic distances vector, and more, in order to prepare the GeoTriMesh object for another round of
     * geodesic computations.
     */
    void reset_distance();
    /**
     * @brief Clears the edge map of the GeoTriMesh.
     * This function clears the edge map of the GeoTriMesh. The edge map is a map of edges to EdgeStructs
     * and is used during geodesic computation.
     */
    void clear_edge_map();

    /**
     * @brief Computes the geodesic distance and path from a selected vertex to all other vertices.
     * This function computes the geodesic distance and path from a selected vertex to all other vertices
     * using the Dijkstra's algorithm. The computed distance and path are returned in the given output
     * arguments `distance` and `path`, respectively.
     * @param selected_v The index of the selected vertex.
     * @param distance The output vector of geodesic distances.
     * @param path The output vector of geodesic paths.
     */
    void compute_geodesic(int selected_v, ScalarVector &distance, KnotVectorVector &path);
    /**
     * @brief Initializes the propagation of geodesic distances from a given vertex.
     * @param selected_v The index of the vertex from which to propagate the geodesic distances.
     * @return true if the initialization was successful, false otherwise.
     */
    bool init_propagation(int selected_v);
    /**
     * @brief Propagates one interval from the priority queue to update the geodesic distance field and path.
     * @return true if an interval was successfully propagated, false if the queue is empty.
     */
    bool propagate_once();
    /**
     * @brief Applies the last step of interval propagation.
     * This function applies the last step of interval propagation that was stored in `last_step`.
     * The modified intervals are updated in the `edge_map`, the interval heap is updated, and
     * the distance field is updated as well.
     */
    void apply_last_step();
    /**
     * @brief Traces the geodesic paths on the mesh.
     * This function traces the geodesic paths on the mesh and stores the knots in `geo_path`.
     */
    void backtracing();

    /**
     * @brief Creates an interval for a given edge and a distance from a source point.
     *
     * @param e The edge for which the interval is created.
     * @param sigma The value of sigma for the interval.
     * @param source The source point.
     * @param e_l The length of the edge. Defaults to -1 to compute the length of the edge.
     * @return The created interval.
     */
    Interval *create_init_interval(Handle e, Vec3 &source);
    /**
     * @brief Creates an interval for a given edge and two parameter values with respect to the edge.
     *
     * @param e The edge for which the interval is created.
     * @param sigma The value of sigma for the interval.
     * @param b0 The lower bound for the parameter value.
     * @param b1 The upper bound for the parameter value.
     * @param source The source point.
     * @param inverted True if the interval is inverted.
     * @param e_l The length of the edge. Defaults to -1 to compute the length of the edge.
     * @return The created interval.
     */
    Interval *create_interval(Handle e, double sigma, Vec3 &source, double e_l = -1);
    /**
     * @brief Creates an interval for a given edge, two parameter values, two distances, and a value of sigma.
     *
     * @param e The edge for which the interval is created.
     * @param b0 The lower bound for the parameter value.
     * @param b1 The upper bound for the parameter value.
     * @param d0 The distance between the source point and the edge point at b0.
     * @param d1 The distance between the source point and the edge point at b1.
     * @param sigma The value of sigma for the interval.
     * @param inverted True if the interval is inverted.
     * @return The created interval.
     */
    Interval *create_interval(Handle e, double sigma, double b0, double b1, Vec3 &source, bool inverted, double e_l = -1);
    /**
     * Creates a new interval object for the given edge with the specified parameters.
     * @param e The edge to create the interval for.
     * @param b0 The left boundary of the interval.
     * @param b1 The right boundary of the interval.
     * @param d0 The distance from the source point to the left boundary of the interval.
     * @param d1 The distance from the source point to the right boundary of the interval.
     * @param sigma The current value of the function sigma.
     * @param inverted Flag indicating whether the edge is inverted.
     * @return A pointer to the newly created interval object.
     */
    Interval *create_interval(Handle e, double b0, double b1, double d0, double d1, double sigma, bool inverted = false);
    /**
     * @brief Calculates the distance value from the source point at the given interval position
     *
     * @param iv Pointer to the interval to calculate the distance value on
     * @param inter Interval position to calculate the distance value at
     * @return double Distance value from the source point at the given interval position
     */
    void update_heap(Interval *from_iv, Propagation &propag);
    /**
     * @brief Updates the intervals in the interval queue heap based on the given altered intervals list
     *
     * @param from_iv Pointer to the interval that the altered intervals list originated from
     * @param altered The list of altered intervals to apply on the intervals in the heap
     */
    void update_heap(Interval *from_iv, IntervalList &altered);
    /**
     * Update the scalar field with the distances of the interval ends.
     *
     * @param to_e The edge to update the scalar field on.
     * @param ivs The set of intervals on the edge.
     */
    void update_sfield(Handle to_e, EdgeStruct &ivs);

    /**
     * Propagates an interval along an edge and updates the interval structure of the mesh.
     * @param iv Interval to propagate
     * @param new_ivs Map to store new intervals created during propagation
     */
    void propagate_interval(Interval *iv, EdgeStructMap &new_ivs);
    /**
     * Merges an interval with a list of old intervals.
     * @param new_iv The new interval to be merged.
     * @param old_ivs A list of old intervals to be merged with the new interval.
     * @param propag A propagation object for tracking changes to the intervals.
     */
    void merge_interval(Interval *new_iv, EdgeStruct &old_ivs, Propagation &propag);

    /**
     * Propagates an interval using case 1.
     * Case 1 is the most simple propagation case. It occurs when the interval intersects the
     * border at two different edges on the same face, and the intersection on the right edge.
     * (b1_it) is after the intersection on the left edge (b0_it).
     * @param b0_it The intersection structure for the left edge.
     * @param b1_it The intersection structure for the right edge.
     * @param iv The interval to be propagated.
     * @param new_ivs The map of newly created intervals.
     */
    void propagate_case1(const InterStruct& b0_it, const InterStruct& b1_it, Interval* iv,
                                EdgeStructMap& new_ivs);
    /**
     * Propagates an interval using case 2.
     * Case 2 occurs when the interval intersects the border at two different edges on the same face,
     * but the two intersections occur on the same edge. In this case, two new intervals are created:
     * one on the left side of the first intersection and another on the right side of the second intersection.
     * @param b0_it The intersection structure for the left edge.
     * @param b1_it The intersection structure for the right edge.
     * @param mid_dist The distance between the left and right edges.
     * @param iv The interval to be propagated.
     * @param new_ivs The map of newly created intervals.
     */
    void propagate_case2(const InterStruct& b0_it, const InterStruct& b1_it, double mid_dist,
                                Interval* iv, EdgeStructMap& new_ivs);
    /**
     * Propagate an interval from one face to an adjacent face through an edge. Handles the case where the interval intersects
     * the border of both faces at the start of the first face. This case produces two new intervals with a pseudo source
     * at the start of the previous edge.
     * @param b0_it The intersection with the border on the first face.
     * @param b1_it The intersection with the border on the second face.
     * @param from_e The edge from which the interval is propagated.
     * @param s The source point of the interval.
     * @param iv The interval to propagate.
     * @param new_ivs Map to store the new intervals created during the propagation.
     */
    void propagate_case3_left(const InterStruct& b0_it, const InterStruct& b1_it, Handle from_e,
                                const Vec2& s, Interval* iv, EdgeStructMap& new_ivs);
    /**
     * @brief Propagates an interval according to case 3 on the right side.
     * This function propagates an interval in case the right endpoint reaches the end of the given edge,
     * and there is a gap between the end of that edge and the start of the next one.
     * The function creates an interval that goes from the given edge's endpoint to the start of the next edge,
     * and two more intervals to cover the gap.
     * @param b0_it The intersecting point at the left side of the edge.
     * @param b1_it The intersecting point at the right side of the edge.
     * @param from_e The edge that the interval is propagated from.
     * @param s The point on the edge from which the interval is propagated.
     * @param iv The interval to be propagated.
     * @param new_ivs A map to store newly created intervals.
     */
    void propagate_case3_right(const InterStruct& b0_it, const InterStruct& b1_it, Handle from_e,
                                const Vec2& s, Interval* iv, EdgeStructMap& new_ivs);

    /**
     * Check if an interval ends at the start of an edge.
     *
     * @param tau Boolean value indicating whether to use the b0 or b1 value of the interval.
     * @param iv The interval to check.
     * @param e_length The length of the edge.
     * @return True if the interval ends at the start of the edge, false otherwise.
     */
    inline bool is_iv_b0_end(bool tau, Interval *iv, double e_length);
    /**
     * Check if an interval ends at the end of an edge.
     *
     * @param tau Boolean value indicating whether to use the b0 or b1 value of the interval.
     * @param iv The interval to check.
     * @param e_length The length of the edge.
     * @return True if the interval ends at the end of the edge, false otherwise.
     */
    inline bool is_iv_b1_end(bool tau, Interval *iv, double e_length);
    /**
     * @brief Returns the s-coordinates of the given interval.
     *
     * Computes the s-coordinates of the given interval and returns it as a Vec2.
     *
     * @param iv Interval to get the s-coordinates from.
     * @return Vec2 containing the s-coordinates of the interval.
     */
    Vec2 get_s_coord(const Interval *iv) const;
    /**
     * @brief Returns the s-coordinates of the given interval using straight line distance.
     *
     * Computes the s-coordinates of the given interval using straight line distance and
     * returns it as a Vec2.
     *
     * @param iv Interval to get the s-coordinates from.
     * @return Vec2 containing the s-coordinates of the interval.
     */
    Vec2 get_s_coord_straight(const Interval *iv) const;
    /**
     * @brief Fills a gap in the diamond shape created by two adjacent edges.
     * This function creates a new interval to fill the gap between the two edges of a diamond shape.
     * @param iv The interval to be propagated.
     * @param cur_e The edge of the diamond shape that the gap is filled in.
     * @param source The source vertex of the propagation.
     * @param cur_sigma The current value of sigma.
     * @param e_length The length of the current edge.
     * @param new_ivs A map to store newly created intervals.
     * @param start A boolean value indicating whether this is the first gap in the diamond shape.
     */    
    void cover_diamond_shape(Interval *iv, Handle cur_e, Vec3 &source, double cur_sigma, double e_length,
                             EdgeStructMap &new_ivs, bool start);
    /**
     * Fills the gaps between two intervals of a diamond-shaped set of edges with new intervals.
     * @param cur_e Handle to the current edge to fill the gaps.
     * @param cur_sigma Double value representing the current sigma value.
     * @param source Vec3 object representing the source point of the current interval.
     * @param new_ivs EdgeStructMap object representing a mapping between edges and intervals.
     * The new intervals will be added to this mapping.
     */
    void fill_diamond_gap(Handle cur_e, double cur_sigma, Vec3 &source, EdgeStructMap &new_ivs);
    /**
     * @brief Projects a 3D point to a 2D point on a face.
     *
     * @param e        An edge on the face.
     * @param to_face  The index of the face.
     * @param p        The 3D point to project.
     * @return         The 2D projection of the point on the face.
     */
    Vec2 project_p(Handle e, int to_face, const Vec3 &p) const;
    /**
     * @brief Unprojects a 2D point on a face to a 3D point.
     *
     * @param e      An edge on the face.
     * @param to_face  The index of the face.
     * @param s_xy     The 2D point to unproject.
     * @return         The 3D point in the unprojected space.
     */
    Vec3 unproject_p(Handle e, int to_face, const Vec2 &s) const;
    /**
     * Tests if a given ray intersects a border edge of a given face of the mesh.
     * @param face The face index to test for intersection.
     * @param from_e The edge from which the search should start.
     * @param coord_map A map from vertex indices to their positions in 2D space.
     * @param ray The ray to test for intersection.
     * @param inter The InterStruct instance that will be updated with information about the intersection if one is found.
     * @return True if an intersection is found, false otherwise.
     */
    bool intersect_face_border(int face, Handle from_e, map<int, Vec2> &coorc_map, Ray<Vec2> ray, InterStruct &inter);
    /**
     * Intersects the border of a given face with a given ray, allowing for a loose intersection (i.e., intersection
     * within a certain distance of the border).
     * @param face The ID of the face to intersect.
     * @param from_e The edge to start from.
     * @param coord_map A map of vertex IDs to their 2D coordinates.
     * @param ray The ray to intersect the face border with.
     * @param inter The resulting intersection structure.
     * @return true if the intersection was successful, false otherwise.
     */    
    bool intersect_face_border_loose(int face, Handle from_e, map<int, Vec2> &coord_map, Ray<Vec2> ray, InterStruct &inter);
    /**
     * @brief Check if two points are close enough (Euclidean distance).
     * @param p0 First point.
     * @param p1 Second point.
     * @return True if the points are close enough, false otherwise.
     */    
    inline bool close_enough(const Vec2 &p0, const Vec2 &p1);
    /**
     * @brief Check if a point is close enough to a ray (distance to the line).
     * @param p The point.
     * @param ray The ray.
     * @return True if the point is close enough to the ray, false otherwise.
     */    
    inline bool close_enough(const Vec2 &p, const Ray<Vec2> &ray);

    /**
     * @brief Get the range of intersection between two intervals.
     * @param v0 The first interval to check.
     * @param v1 The second interval to check.
     * @return Range The range of intersection.
     */
    Range get_intersection_range(const Interval *v0, const Interval *v1) const;
    /**
     * @brief Compare two intervals' distances within a given range.
     * @param v0 The first interval to compare.
     * @param v1 The second interval to compare.
     * @param range The range to compare the intervals' distances.
     * @return short -1 if v0 has shorter distance all over the range, 1 if v1 has shorter distance all over the range, 0 otherwise.
     */
    short absolute_compare(const Interval *v0, const Interval *v1, const Range &range) const;
    /**
     * @brief Merge two intervals where the first interval is shorter than the second interval
     * @param v0 The first interval
     * @param v1 The second interval
     * @param propag The propagation object
     */
    void absolute_merging_v0_shorter_v1(Interval *v0, Interval *v1, Propagation &propag);
    /**
     * @brief Merge two intervals where the second interval is shorter than the first interval
     * @param v0 The first interval
     * @param v1 The second interval
     * @param propag The propagation object
     * @param rg The range of the overlap between the two intervals
     */
    void absolute_merging_v1_shorter_v0(Interval *v0, Interval *v1, Propagation &propag, Range &rg);
    /**
     * @brief Computes the x-coordinate of the intersection between two intervals on the x-axis
     * @param v0 The first interval
     * @param v1 The second interval
     * @param rg The range of the overlap between the two intervals
     * @return The x-coordinate of the intersection point
     */
    double get_px_value(Interval *v0, Interval *v1, Range &rg);
    /**
     * @brief Merge two intervals partially if they overlap in a range.
     * @param v0 Interval to be merged.
     * @param v1 Interval to be merged.
     * @param propag Propagation object to keep track of changes.
     * @param rg Range where the overlap occurs.
     */
    void partial_merging(Interval *v0, Interval *v1, Propagation &propag, Range &rg);
    /**
     * @brief Tests if v0->get_b0() <= range.first, and returns 4 if true, 0 otherwise.
     * @param v0 Pointer to the first Interval.
     * @param v1 Pointer to the second Interval.
     * @param range Reference to the range object.
     * @return An unsigned char value. Returns 4 if v0->get_b0() <= range.first, 0 otherwise.
     */
    unsigned char test001(Interval *v0, Interval *v1, Range &range);
    /**
     * @brief Tests if range.second < v0->get_b1(), and returns 2 if true, 0 otherwise.
     * @param v0 Pointer to the first Interval.
     * @param v1 Pointer to the second Interval.
     * @param range Reference to the range object.
     * @return An unsigned char value. Returns 2 if range.second < v0->get_b1(), 0 otherwise.
     */
    unsigned char test010(Interval *v0, Interval *v1, Range &range);
    /**
     * @brief Tests if distance_at(v0, range.first) < distance_at(v1, range.first), and returns 1 if true, 0 otherwise.
     * @param v0 Pointer to the first Interval.
     * @param v1 Pointer to the second Interval.
     * @param range Reference to the range object.
     * @return An unsigned char value. Returns 1 if distance_at(v0, range.first) < distance_at(v1, range.first), 0 otherwise.
     */
    unsigned char test100(Interval *v0, Interval *v1, Range &range);
    /**
     * @brief Performs post-processing on a Propagation object.
     * @param propag Reference to the Propagation object to be post-processed.
     * @param b Pointer to the Interval object.
     */
    void post_processing(Propagation &propag, Interval *b);

    // auxiliarly geometric functions
    /**
     * @brief Returns the length of an edge.
     * @param e Handle to the edge.
     * @return double Length of the edge.
     */
    double edge_length(const Handle &e) const;
    /**
     * @brief Returns a unit vector representing the direction of an edge.
     * @param e Handle to the edge.
     * @return Vec3 Unit vector representing the direction of the edge.
     */
    Vec3 edge_vector(const Handle &e) const;
    /**
     * @brief Computes the 3D point at a given interval along an edge.
     * @param handle Handle to the edge.
     * @param interv Interval along the edge, as a ratio of its length.
     * @return Vec3 The 3D point at the given interval along the edge.
     */
    Vec3 edge_point(const Handle &handle, double interv) const;
    /**
     * @brief Computes the 3D point at a given ratio along an edge.
     * @param handle Handle to the edge.
     * @param ratio Ratio along the edge.
     * @return Vec3 The 3D point at the given ratio along the edge.
     */
    Vec3 edge_point_unit(const Handle &handle, double ratio) const;
    /**
     * @brief Computes the interval on an edge that is closest to a given 2D point.
     * @param p0 First end point of the edge.
     * @param p1 Second end point of the edge.
     * @param point The 2D point.
     * @return double The interval on the edge that is closest to the point.
     */
    double edge_interval(const Vec2 &p0, const Vec2 &p1, const Vec2 &point) const;
    /**
     * @brief Computes the ratio on an edge that is closest to a given 2D point.
     * @param p0 First end point of the edge.
     * @param p1 Second end point of the edge.
     * @param point The 2D point.
     * @return double The ratio on the edge that is closest to the point.
     */
    double edge_interval_unit(const Vec2 &p0, const Vec2 &p1, const Vec2 &point) const;
    /**
     * @brief Calculates the distance from the source point at the given interval
     *
     * @param iv Pointer to the interval to calculate the distance on
     * @param inter Interval position to calculate the distance at
     * @return double Distance from the source point at the given interval
     */
    double distance_at(const Interval *iv, const double inter) const;
    /**
     * @brief Calculates the distance value from the source point at the given interval position
     *
     * @param iv Pointer to the interval to calculate the distance value on
     * @param inter Interval position to calculate the distance value at
     * @return double Distance value from the source point at the given interval position
     */
    double dvalue_at(const Interval *iv, const double inter) const;

    /**
     * @brief Sets the error threshold value for the GeoPath algorithm
     * @param thres The error threshold value to set
     */
    void set_error_threshold(double thres);

    /**
     * @brief Computes the normal of a vertex by averaging the normals of the faces it belongs to
     * @param vid The index of the vertex to compute the normal for
     * @return Vec3 The normal of the vertex
     */
    Vec3 vertex_normal(int vid);
};

extern bool DEBUG;

#endif