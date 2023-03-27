#include "trimesh.hpp"
#include "geom3d.hpp"

#include <fstream>

void TriMesh::clear()
{
    m_vertex.clear();
    m_face.clear();
    m_fnormal.clear();
}

void TriMesh::initialize()
{
    m_selected_vertex = -1;
    m_selected_face = -1;

    compute_fnormal();
//    initVtxFaceMap();
}

void TriMesh::compute_bbox(Vec3& min, Vec3& max)
{
    if( !m_vertex.size() )  min = max = 0.0;
    else                    min = max = m_vertex[0];

    for(int i=1; i<m_vertex.size(); i++)
    for(int j=0; j<3; j++)
    {
        if( m_vertex[i][j] < min[j])  min[j]=m_vertex[i][j];
        if( m_vertex[i][j] > max[j])  max[j]=m_vertex[i][j];
    }
}

void TriMesh::compute_fnormal()
{
    m_fnormal.clear();
    for (FaceList::iterator iter = m_face.begin(); iter != m_face.end(); iter++) 
    {
        const Face& f = *iter;
        m_fnormal.push_back(triangle_normal(
                    m_vertex[f[0]], m_vertex[f[1]], m_vertex[f[2]]));
    }
}

void TriMesh::normalize()
{
    Vec3 bmin, bmax;
    compute_bbox(bmin, bmax);
    double d = norm(bmax-bmin);
    Vec3 center = bmin + (bmax - bmin)/2.;

    int vsize = m_vertex.size();
    for (int i=0; i<vsize; i++)
        m_vertex[i] = 1/d * (m_vertex[i] - center);

    compute_fnormal();
}

void TriMesh::read_from_file(const string& filename)
{
    ifstream fin(filename.c_str());
    if (!fin) {
        cerr << "Cannot open file " << filename << endl;
        return;
    }

    this->clear();
    char line[1024];
    while (fin.getline(line, 1024)) {
        if (line[0] == 'v') {
            Vec3 v;
            sscanf(line, "v %lf %lf %lf", &v[0], &v[1], &v[2]);
            m_vertex.push_back(v);
        }
        else if (line[0] == 'f') {
            Face f;
            sscanf(line, "f %d %d %d", &f[0], &f[1], &f[2]);
            m_face.push_back(f);
        }
    }
    this->compute_fnormal();

    fin.close();
}

void TriMesh::write_to_file(const string& filename)
{
    ofstream fout(filename.c_str());
    if (!fout) {
        cerr << "Cannot open file " << filename << endl;
        return;
    }

    for (int i = 0; i < m_vertex.size(); i++) {
        fout << "v " << m_vertex[i][0] << " " << m_vertex[i][1] << " " << m_vertex[i][2] << endl;
    }
    for (int i = 0; i < m_face.size(); i++) {
        fout << "f " << m_face[i][0] << " " << m_face[i][1] << " " << m_face[i][2] << endl;
    }

    fout.close();
}