#include "mls.hpp"
#include "auxiliar.hpp"
#include "general_math.hpp"

MLS::MLS(GeoTriMesh *m)
: mesh(m)
{
	int vsize = mesh->m_vertex.size();
}

MLS::~MLS()
{
    if (mesh) delete mesh;
	clear_distances();
}

void MLS::clear_distances()
{
    LOG(INFO) << "MLS::clear_distances()";
    for (auto& path : paths)
        for (auto& knot : path)
            knot.clear();
    paths.clear(); 
    for (auto& dist : distances)
        dist.clear();
	distances.clear();
}

void MLS::compute_distances(int selected_v)
{
    LOG(INFO) << "MLS::compute_distances(" << selected_v << ")";
	clear_distances();

	int vsize = mesh->m_vertex.size();

    // geodesic distance
    ScalarVector dist;
    KnotVectorVector path;
    mesh->compute_geodesic(selected_v, dist, path);
    distances.push_back(dist);
    paths.push_back(path);

    cout << "finished computing geodesic for vertex " << selected_v << endl;
}


double MLS::project_distance(KnotVector& kv, VectorVector& field)
{
    LOG(INFO) << "MLS::project_distance()";
	int ksize = kv.size(), index=1;

	if (ksize < 2)
		return 0.0;

	double acc_dist=0, proj_dist=0;
	KnotPair knext, kprev = kv[0];
	Vec3 vnext, vprev = mesh->edge_point(kprev.first, kprev.second);
	Vec3 vec, grad;

	do {
		knext = kv[index];
		vnext = mesh->edge_point(knext.first, knext.second);
		vec = vnext-vprev;
		grad = interpolate_gradient(field, (Handle) kprev.first, kprev.second);
		proj_dist = vec * grad;
		acc_dist += proj_dist;

		index++;
		kprev = knext;
		vprev = vnext;
	} 
	while (index < ksize);

	if (acc_dist < 0.0)
		acc_dist = -acc_dist;

	return acc_dist;
}

Vec3 MLS::interpolate_gradient(VectorVector& field, Handle e, double interval)
{
    LOG(INFO) << "MLS::interpolate_gradient()";
	double alpha = interval / mesh->edge_length(e);
	Vec3 inter = (1.0-alpha) * field[e->Org()] + alpha * field[e->Dest()];

	// TODO:: this could be removed
	unitize(inter);

	return inter;
}

void MLS::print_knot_vector(const KnotVector& kn)
{
	KnotVector::const_iterator kit = kn.begin();
	for ( ; kit != kn.end(); kit++) 
		cout << "k: " << *(*kit).first << "  " << (*kit).second << endl;
}
