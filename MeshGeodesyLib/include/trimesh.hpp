#ifndef TRIMESH_INCLUDED
#define TRIMESH_INCLUDED

#include "vec3.hpp"

#include <vector>
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

public:
    VertexList vertex;
    FaceList face;
    NormalVector fnormal;

  	int selected_vertex;
	int selected_face;

public:
    TriMesh() {}
    TriMesh(const TriMesh& m) { *this = m; }
    TriMesh(VertexList& vlist, FaceList& flist) 
        : vertex(vlist), face(flist) {}
    ~TriMesh() { clear(); }
    TriMesh& operator=(const TriMesh& m)
    {
        vertex = m.vertex;
        face = m.face;
        fnormal = m.fnormal;
        return *this;
    }
    void clear();
	void initialize();

    void computeBBox(Vec3& min, Vec3& max);
	void computeFNormal();
	void normalize();

	// in/output methods
	void read_from_file(const string& filename);
	void write_to_file(const string& filename);
};
// TRIMESH_INCLUDED
#endif
