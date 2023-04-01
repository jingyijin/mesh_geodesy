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

void MLS::print_knot_vector(const KnotVector& kn)
{
    for (auto& k : kn)
        cout << "k: " << *k.first << "  " << k.second << endl;
}
