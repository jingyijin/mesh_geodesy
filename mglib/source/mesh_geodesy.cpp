/************************************************************************
 * File description: Mesh Geodesy class as the main class to store distances
 * and paths for the geodesic computation.
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "mesh_geodesy.hpp"
#include "general_math.hpp"

#include <fstream>
#include <algorithm>
#include <iomanip>

MeshGeodesy::MeshGeodesy(GeoTriMesh *m)
: m_mesh(m), m_source_v(-1)
{
    int vsize = m_mesh->m_vertex.size();
}

MeshGeodesy::~MeshGeodesy()
{
    if (m_mesh) delete m_mesh;
    clear_distances();
}

void MeshGeodesy::clear_distances()
{
    LOG(INFO) << "MeshGeodesy::clear_distances()";
    for (auto& path : m_path_points)
        path.clear();
    m_path_points.clear(); 
    m_distances.clear();
}

void MeshGeodesy::compute_distances(int selected_v)
{
    LOG(INFO) << "MeshGeodesy::compute_distances(" << selected_v << ")";
    clear_distances();

    int vsize = m_mesh->m_vertex.size();

    m_source_v = selected_v;
    // geodesic distance
    m_mesh->compute_geodesic(selected_v, m_distances, m_path_points);
}

void MeshGeodesy::sort_faces_by_distance()
{
    cout << "MeshGeodesy::sort_faces_by_distance()" << endl;
    // sort faces by distance
    vector<pair<double, Face>> dist_face_pairs;
    int i;
    double min_dist;
    for (auto& f : m_mesh->m_face)
    {
        min_dist = DBL_MAX;
        // min of dist at vertices
        for (i=0; i<3; ++i)
            min_dist = min(min_dist, m_distances[f[i]]);
        min_dist /= 3.;
        dist_face_pairs.push_back(make_pair(min_dist, f));
    }
    sort(dist_face_pairs.begin(), dist_face_pairs.end(), 
        [](const pair<double, Face>& p1, const pair<double, Face>& p2) {
            return p1.first < p2.first;
    });

    // update the face list
    m_mesh->m_face.clear();
    for (auto& p : dist_face_pairs)
        m_mesh->m_face.push_back(p.second);
    m_mesh->compute_fnormal();
}

void MeshGeodesy::save_geodesic(const string& filename)
{
    LOG(INFO) << "MeshGeodesy::save_distances(" << filename << ")";
    ofstream fout(filename);
    if (!fout.is_open()) {
        LOG(ERROR) << "Cannot open file " << filename;
        return;
    }

    // Set precision for floating-point output
    fout << std::fixed << std::setprecision(6);

    // Save the vertices with their distance values
    GeoTriMesh::VertexList& vertex = m_mesh->m_vertex;
    for (size_t i = 0; i < vertex.size(); ++i) {
        fout << "d " << vertex[i][0] << " " << vertex[i][1] << " " << vertex[i][2] << " ";
        fout << (m_distances.empty() ? 0.0 : m_distances[i]) << '\n';
    }

    // Save the faces
    GeoTriMesh::FaceList& face = m_mesh->m_face;
    for (const auto& f : face)
        fout << "f " << f[0] << " " << f[1] << " " << f[2] << '\n';

    // Save the paths
    for (size_t i = 0; i < m_path_points.size(); ++i) {
        const auto& path = m_path_points[i];
        fout << "p " << m_source_v << " " << i << " " << path.size();
        for (const auto& v : path)
            fout << " " << v[0] << " " << v[1] << " " << v[2];
        fout << '\n';
    }

    fout.close();
}

void MeshGeodesy::load_geodesic(const string& filename)
{
    LOG(INFO) << "MeshGeodesy::load_distances(" << filename << ")";
    ifstream fin(filename);
    if (!fin.is_open())
    {
        LOG(ERROR) << "Cannot open file " << filename;
        return;
    }

    m_mesh->clear();
    std::string line;
    GeoTriMesh::VertexList& vertex = m_mesh->m_vertex;
    GeoTriMesh::FaceList& face = m_mesh->m_face;
    m_distances.clear();
    m_path_points.clear();

    while (std::getline(fin, line)) {
        std::istringstream iss(line);
        char type;
        iss >> type;

        if (type == 'd') {
            Vec3 vert;
            double distance;
            iss >> vert[0] >> vert[1] >> vert[2] >> distance;
            vertex.push_back(vert);
            m_distances.push_back(distance);
        }
        else if (type == 'f') {
            Face f;
            iss >> f[0] >> f[1] >> f[2];
            face.push_back(f);
        }
        else if (type == 'p') {
            int source_vertex_id, target_vertex_id, num_points;
            iss >> source_vertex_id >> target_vertex_id >> num_points;
            m_source_v = source_vertex_id;

            GeoTriMesh::VertexList path;
            for (int i = 0; i < num_points; ++i) {
                Vec3 point;
                iss >> point[0] >> point[1] >> point[2];
                path.push_back(point);
            }
            m_path_points.push_back(path);
        }
    }
    fin.close();

    m_mesh->initialize();
    m_mesh->compute_fnormal();
}

