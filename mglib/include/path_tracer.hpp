#ifndef PATH_TRACER_INCLUDED
#define PATH_TRACER_INCLUDED

#include "trimesh.hpp"
#include "interval.hpp"

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
	GeoTriMesh& mesh;
	PathKnotVectorVector& path_knot;

	State_type state;
	int next_vertex, target_vertex, vindex;
	Interval* min_iv;
	Handle min_e, cur_e, prev_e;
	InterStruct inter;
	Vec2 s, start_coord, inter_p2D;
	int to_face;
	Ray<Vec2> ray;
	map<Handle, bool> current_ring;

public:
	PathTracer(GeoTriMesh& m, PathKnotVectorVector& path_knot);

	void clear_path();
	void trace_path();
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