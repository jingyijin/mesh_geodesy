#ifndef MLS_INCLUDED
#define MLS_INCLUDED

#include "geomesh.hpp"

class MLS
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
	vector<ScalarVector> distances;
	vector<KnotVectorVector> paths;

	// constraints
	ValueConstraintMap value_constraints;
	GradientConstraintMap gradient_constraints;
	vector<VectorVector> gradient_field;

	ScalarVector scalar_field;

public:
	MLS(GeoTriMesh *m);
	~MLS();

	void clear();
	void clear_distances();

	void compute_distances(int selected_v);

	// methods related to computing projected path distances
	double project_distance(KnotVector& kv, VectorVector& field);
    Vec3 interpolate_gradient(VectorVector& field, Handle e, double interval);

	// in/output methods
	void print_knot_vector(const KnotVector& kn);
};

#endif