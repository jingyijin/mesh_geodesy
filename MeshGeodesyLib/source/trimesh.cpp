#include "trimesh.hpp"
#include "geom3d.hpp"

#include <fstream>

void TriMesh::clear()
{
	vertex.clear();
	face.clear();
	fnormal.clear();
}

void TriMesh::initialize()
{
	selected_vertex = -1;
	selected_face = -1;

	computeFNormal();
//	initVtxFaceMap();
}

void TriMesh::computeBBox(Vec3& min, Vec3& max)
{
    if( !vertex.size() )  min = max = 0.0;
    else                  min = max = vertex[0];

    for(int i=1; i<vertex.size(); i++)
	for(int j=0; j<3; j++)
	{
	    if( vertex[i][j] < min[j])  min[j]=vertex[i][j];
	    if( vertex[i][j] > max[j])  max[j]=vertex[i][j];
	}
}

void TriMesh::computeFNormal()
{
	fnormal.clear();
	for (FaceList::iterator iter = face.begin(); iter != face.end(); iter++) 
	{
		const Face& f = *iter;
		fnormal.push_back(triangle_normal(vertex[f[0]], vertex[f[1]], vertex[f[2]]));
	}
}

void TriMesh::read_from_file(const string& filename)
{
    ifstream fin(filename.c_str());
    if (!fin) {
        cerr << "Cannot open file " << filename << endl;
        return;
    }

    char line[1024];
    while (fin.getline(line, 1024)) {
        if (line[0] == 'v') {
            Vec3 v;
            sscanf(line, "v %lf %lf %lf", &v[0], &v[1], &v[2]);
            vertex.push_back(v);
        }
        else if (line[0] == 'f') {
            Face f;
            sscanf(line, "f %d %d %d", &f[0], &f[1], &f[2]);
            f -= 1;
            face.push_back(f);
        }
    }
    this->computeFNormal();

    fin.close();
}

void TriMesh::write_to_file(const string& filename)
{
    ofstream fout(filename.c_str());
    if (!fout) {
        cerr << "Cannot open file " << filename << endl;
        return;
    }

    for (int i = 0; i < vertex.size(); i++) {
        fout << "v " << vertex[i][0] << " " << vertex[i][1] << " " << vertex[i][2] << endl;
    }
    for (int i = 0; i < face.size(); i++) {
        fout << "f " << face[i][0] << " " << face[i][1] << " " << face[i][2] << endl;
    }

    fout.close();
}