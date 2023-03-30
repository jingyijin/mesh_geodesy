#ifndef TRIMESH_INCLUDED
#define TRIMESH_INCLUDED

#include "vec3.hpp"

#include <vector>
#include <set>
#include <unordered_map>

using namespace std;
using namespace gfx;

struct Face : public TVec3<id_t>
{
    size_t size() { return 3; }
};

class TriMesh
{
public:
    typedef vector<Vec3> VertexList;
    typedef vector<Face> FaceList;
    typedef vector<Vec3> NormalVector;

	typedef set<id_t> FaceSet;
	typedef set<id_t> VertexSet;

	typedef Face Cell;
	typedef FaceList CellList;

public:
    VertexList m_vertex;
    FaceList m_face;
    NormalVector m_fnormal;

    int m_selected_vertex;
    int m_selected_face;

public:
    TriMesh() {}
    TriMesh(const TriMesh& m) { *this = m; }
    TriMesh(VertexList& vlist, FaceList& flist) 
        : m_vertex(vlist), m_face(flist) {}
    ~TriMesh() { clear(); }
    TriMesh& operator=(const TriMesh& m)
    {
        m_vertex = m.m_vertex;
        m_face = m.m_face;
        m_fnormal = m.m_fnormal;
        return *this;
    }
    void clear();
    void initialize();

    void compute_bbox(Vec3& min, Vec3& max);
    void compute_fnormal();
    void normalize();

	int other_vertex(int fid, int v0, int v1) const;

    // in/output methods
    void read_from_file(const string& filename);
    void write_to_file(const string& filename);
};
// TRIMESH_INCLUDED
#endif
