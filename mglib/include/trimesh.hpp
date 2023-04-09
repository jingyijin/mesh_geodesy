#ifndef TRIMESH_INCLUDED
#define TRIMESH_INCLUDED

/************************************************************************
 * File description: Triangular mesh class
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/


#include "vec3.hpp"

#include <vector>
#include <set>
#include <unordered_map>

using namespace std;
using namespace gfx;

/**
 * @brief A structure representing a triangular face in a mesh.
 */
struct Face : public TVec3<id_t>
{
    /**
     * @brief Returns the size of the triangular face (always 3).
     * @return The number of vertices in the face.
     */
    size_t size() { return 3; }
};

/**
 * @brief A class representing a triangular mesh.
 */
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
    VertexList m_vertex;    /*< List of vertices in the mesh. */
    FaceList m_face;        /*< List of faces in the mesh. */
    NormalVector m_fnormal; /*< List of face normals. */

    int m_selected_vertex;  /*< The index of the selected vertex. */
    int m_selected_face;    /*< The index of the selected face. */

public:
    /**
     * @brief Default constructor.
     */
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
    /**
     * @brief Clears the mesh data.
     */
    void clear();
    /**
     * @brief Initializes the mesh.
     */
    void initialize();

    /**
     * @brief Computes the bounding box of the mesh.
     * @param[out] min The minimum point of the bounding box.
     * @param[out] max The maximum point of the bounding box.
     */
    void compute_bbox(Vec3& min, Vec3& max);
    /**
     * @brief Computes the face normals for the mesh.
     */
    void compute_fnormal();
    /**
     * @brief Normalizes the coordinates of the mesh to be between [0-1].
     */
    void normalize();

    /**
     * @brief Returns the third vertex in the given face.
     * @param fid The face index.
     * @param v0 The index of the first vertex.
     * @param v1 The index of the second vertex.
     * @return The index of the third vertex in the face.
     */
    int other_vertex(int fid, int v0, int v1) const;

    /**
     * @brief Reads the mesh data from a file.
     * @param filename The input file name.
     */
    void read_from_file(const string& filename);
    /**
     * @brief Writes the mesh data to a file.
     * @param filename The output file name.
     */
    void write_to_file(const string& filename);
};
// TRIMESH_INCLUDED
#endif
