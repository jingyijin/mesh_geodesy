/************************************************************************
 * File description: MLS class as the main class to store distances
 * and paths for the geodesic computation.
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "mls.hpp"
#include "general_math.hpp"

#include <fstream>
#include <algorithm>

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
        path.clear();
    paths.clear(); 
    distances.clear();
}

void MLS::compute_distances(int selected_v)
{
    LOG(INFO) << "MLS::compute_distances(" << selected_v << ")";
    clear_distances();

    int vsize = mesh->m_vertex.size();

    // geodesic distance
    mesh->compute_geodesic(selected_v, distances, paths);

    cout << "finished computing geodesic for vertex " << selected_v << endl;
}

void MLS::sort_faces_by_distance()
{
    cout << "MLS::sort_faces_by_distance()" << endl;
    // sort faces by distance
    vector<pair<double, Face>> dist_face_pairs;
    int i;
    double min_dist;
    for (auto& f : mesh->m_face)
    {
        min_dist = DBL_MAX;
        // min of dist at vertices
        for (i=0; i<3; ++i)
            min_dist = min(min_dist, distances[f[i]]);
        min_dist /= 3.;
        dist_face_pairs.push_back(make_pair(min_dist, f));
    }
    sort(dist_face_pairs.begin(), dist_face_pairs.end(), 
        [](const pair<double, Face>& p1, const pair<double, Face>& p2) {
            return p1.first < p2.first;
    });

    // update the face list
    mesh->m_face.clear();
    for (auto& p : dist_face_pairs)
        mesh->m_face.push_back(p.second);
    mesh->compute_fnormal();
}

void MLS::save_distances(const string& filename)
{
    LOG(INFO) << "MLS::save_distances(" << filename << ")";
    ofstream fout(filename);
    if (!fout.is_open())
    {
        LOG(ERROR) << "Cannot open file " << filename;
        return;
    }

    for (auto& dist : distances)
        fout << "s " << dist << endl;

    fout.close();
}

void MLS::load_distances(const string& filename)
{
    LOG(INFO) << "MLS::load_distances(" << filename << ")";
    ifstream fin(filename);
    if (!fin.is_open())
    {
        LOG(ERROR) << "Cannot open file " << filename;
        return;
    }

    distances.clear();
    // load distances
    char line[1024];
    double dist=0.;
    while (fin.getline(line, 1024)) {
        if (line[0] == 's') {
            sscanf(line, "s %lf", &dist);
            distances.push_back(dist);
        }
    }
    fin.close();
}

void MLS::print_knot_vector(const KnotVector& kn)
{
    for (auto& k : kn)
        cout << "k: " << *k.first << "  " << k.second << endl;
}
