#ifndef GEO_SMF_INCLUDED
#define GEO_SMF_INCLUDED

#include "trimesh.hpp"
#include "edgestruct.hpp"
#include "manifold.hpp"
#include "heap.hpp"
#include "ray.hpp"
#include "propagation.hpp"
#include "typedef.hpp"
#include "last_step.hpp"

#include <map>

typedef ManifoldGraph<TriMesh> ManifoldGraphT;
typedef ManifoldGraphT::Handle Handle;

class GeoTriMesh : public TriMesh
{
public:
    typedef set<Handle> HandleSet;
    typedef EdgeStruct::RangeSet RangeSet;
    typedef vector<double> ScalarVector;
    typedef pair<int, int> IntPair;

    struct cmp_handle
    {
        bool operator()(const Handle a, const Handle b) const
        {
            if (a->Org() == b->Org())
                return a->Dest() < b->Dest();
            return a->Org() < b->Org();
        }
    };
    struct cmp_intpair
    {
        bool operator()(const IntPair& a, const IntPair& b) const
        {
            if (a.first == b.first)
                return a.second < b.second;
            return a.first < b.first;
        }
    };
    struct InterStruct
    { 
        Handle e; 
        double inter;
        double t;

        InterStruct() { e = NULL; inter = 0.; t = 0.; }
        bool is_interp_end(double e_value) const { return FEQ(inter, e_value); }
    };

    typedef map<Handle, EdgeStruct, cmp_handle> EdgeStructMap;
    typedef map<IntPair, RangeSet, cmp_intpair> EdgeCoveringMap;

    typedef pair<Handle, Interval*> IntervalPair;

    typedef pair<Handle, double> PathKnotPair;
    typedef vector<PathKnotPair> PathKnotVector;
    typedef vector<PathKnotVector> KnotVectorVector;

public:
    ManifoldGraphT* graph;

    // ******************************
    // for geodesic computation
    MxHeap i_queue;
    EdgeStructMap edge_map;
    ScalarVector geo_distance;
    KnotVectorVector geo_path;

    int prop_step_size;
    int num_prop_edge;
    int stop_prop_number;
    bool is_first_propagation;
    double ERROR_TOLERANCE;

    LastStepInfo last_step;

    // for checking purpose
    EdgeCoveringMap edge_covering;

public:
    GeoTriMesh();
    GeoTriMesh(TriMesh *m);
    ~GeoTriMesh();

    void initialize();
    void reset_distance();
    void clear_edge_map();

    // core computation
    void compute_geodesic(int selected_v, ScalarVector& distance, KnotVectorVector& path);
    bool init_propagation(int selected_v);
    bool propagate_once();
    void apply_last_step();
    void backtracing();

    Interval* create_init_interval(Handle e, Vec3& source);
    Interval* create_interval(Handle e, double sigma, Vec3& source, double e_l=-1);
    Interval* create_interval(Handle e, double sigma, double b0, double b1, Vec3& source, bool inverted, double e_l=-1);
    Interval* create_interval(Handle e, double b0, double b1, double d0, double d1, double sigma, bool inverted=false);
    void update_heap(Interval* from_iv, Propagation& propag);
    void update_heap(Interval* from_iv, IntervalList& altered);
    void update_sfield(Handle to_e, EdgeStruct& ivs);

    // propagation and merging of the intervals
    void propagate_interval(Interval* iv, EdgeStructMap& new_ivs);
    void merge_interval(Interval* new_iv, EdgeStruct& old_ivs, Propagation& propag);

    // propagation cases
    void propagate_case1(const InterStruct& b0_it, const InterStruct& b1_it, Interval* iv, 
                                EdgeStructMap& new_ivs);
    void propagate_case2(const InterStruct& b0_it, const InterStruct& b1_it, double mid_dist,
                                Interval* iv, EdgeStructMap& new_ivs);
    void propagate_case3_left(const InterStruct& b0_it, const InterStruct& b1_it, Handle from_e, 
                                const Vec2& s, Interval* iv, EdgeStructMap& new_ivs);
    void propagate_case3_right(const InterStruct& b0_it, const InterStruct& b1_it, Handle from_e, 
                                const Vec2& s, Interval* iv, EdgeStructMap& new_ivs);

    // propagation related methods
    bool is_iv_b0_end(bool tau, Interval* iv, double e_length);
    bool is_iv_b1_end(bool tau, Interval* iv, double e_length);
    Vec2 get_s_coord(const Interval* iv) const;
    Vec2 get_s_coord_straight(const Interval* iv) const;
    void cover_diamond_shape(Interval* iv, Handle cur_e, Vec3& source, double cur_sigma, double e_length, 
                                EdgeStructMap& new_ivs, bool start);
    void fill_diamond_gap(Handle cur_e, double cur_sigma, Vec3& source, EdgeStructMap& new_ivs);
    Vec2 project_p  (Handle e, int to_face, const Vec3& p) const;
    Vec3 unproject_p(Handle e, int to_face, const Vec2& s) const;
    bool intersect_face_border(int face, Handle from_e, map<int,Vec2>& coorc_map, Ray<Vec2> ray, InterStruct& inter);
    bool intersect_face_border_loose(int face, Handle from_e, map<int,Vec2>& coord_map, Ray<Vec2> ray, InterStruct& inter);
    bool close_enough(const Vec2& p0, const Vec2& p1);
    bool close_enough(const Vec2& p, const Ray<Vec2>& ray);

    // merging related methods
    Range get_intersection_range(const Interval* v0, const Interval* v1) const;
    short absolute_compare(const Interval* v0, const Interval* v1, const Range& range) const;
    void absolute_merging_v0_shorter_v1(Interval* v0, Interval* v1, Propagation& propag);
    void absolute_merging_v1_shorter_v0(Interval* v0, Interval* v1, Propagation& propag, Range& rg);
    double get_px_value(Interval* v0, Interval* v1, Range& rg);
    void partial_merging(Interval* v0, Interval* v1, Propagation& propag, Range& rg);
    unsigned char test001(Interval* v0, Interval* v1, Range& range);
    unsigned char test010(Interval* v0, Interval* v1, Range& range);
    unsigned char test100(Interval* v0, Interval* v1, Range& range);
    void post_processing(Propagation& propag, Interval* b);

    // approximation methods
    void approximate_interval(Propagation& propag);
    Vec3 create_line(const Vec2& p0, const Vec2& p1) const;
    bool is_mergeble(Interval* a, Interval* b, const Vec2& s0, const Vec2& s1, Vec2& new_s, double& new_sigma);
    double error_tolerance(Interval* a, Interval* b, const Vec2& new_s, const double new_sigma);
    void get_dp_interp(double ext0, double ext1, double A, double B, double C, vector<Vec2>& interp);
    void get_dv_interp(const Vec3& L0, const Vec3& L1, double A, double B, double C, vector<Vec2>& dv_interp);
    pair<Vec2, double> min_sigma(vector<Vec2>& interp, double alpha, double beta);
    bool test_visibility(const vector<Vec2>& interp, Vec3& l0, Vec3& l1, double b0, double b1, vector<Vec2>& res_interp);
    bool test_positivity(const vector<Vec2>& interp, double ext0, double ext1, vector<Vec2>& res_interp);

    // auxiliarly geometric functions
    double edge_length(const Handle& e) const;
    Vec3 edge_vector(const Handle& e) const;
    Vec3 edge_point(const Handle& handle, double interv) const;
    Vec3 edge_point_unit(const Handle& handle, double ratio) const;
    double edge_interval(const Vec2& p0, const Vec2& p1, const Vec2& point) const;
    double edge_interval_unit(const Vec2& p0, const Vec2& p1, const Vec2& point) const;
    double distance_at(const Interval* iv, const double inter) const;
    double dvalue_at(const Interval* iv, const double inter) const;

    // debug methods
    void check_edge_covering();
    bool check_validity(EdgeStruct& ivs);
    bool check_validity(IntervalList& ivs);
    bool has_invalid_s_coord(const IntervalList& es) const;

    void set_error_threshold(double thres);

    // general operations
    Vec3 vertex_normal(int vid);
};

extern bool DEBUG;

#endif