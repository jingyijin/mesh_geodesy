#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <string>
#include <glog/logging.h>

#include "geomesh.hpp"
#include "mesh_geodesy.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    google::InitGoogleLogging(argv[0]);

    int opt;
    string input_file, output_file;
    int source_vertex = -1, target_vertex = -1;

    static struct option long_options[] = {
        {"input_file", required_argument, nullptr, 'i'},
        {"source_vertex", required_argument, nullptr, 's'},
        {"target_vertex", required_argument, nullptr, 't'},
        {"output_file", required_argument, nullptr, 'o'},
        {nullptr, 0, nullptr, 0}
    };

    while ((opt = getopt_long(argc, argv, "i:s:t:o:", long_options, nullptr)) != -1)
    {
        switch (opt)
        {
        case 'i':
            input_file = optarg;
            break;
        case 's':
            source_vertex = std::stoi(optarg);
            break;
        case 't':
            target_vertex = std::stoi(optarg);
            break;
        case 'o':
            output_file = optarg;
            break;
        default:
            std::cerr << "Usage: " << argv[0] << " [-i|--input_file <filename>] [-s|--source_vertex <vertex_id>] "
                      << "[-t|--target_vertex <vertex_id>] [-o|--output_file <filename>]" << std::endl;
            return EXIT_FAILURE;
        }
    }

    // Check if the required arguments are provided and valid
    if (input_file.empty() || output_file.empty() || source_vertex < 0)
    {
        std::cerr << "Usage: " << argv[0] << " [-i|--input_file <filename>] [-s|--source_vertex <vertex_id>] "
                  << "[-t|--target_vertex <vertex_id>] [-o|--output_file <filename>]" << std::endl;
        return EXIT_FAILURE;
    }

    try
    {
        TriMesh *tri = new TriMesh();
        tri->read_from_file(input_file);
        tri->initialize();
        GeoTriMesh *m_mesh = new GeoTriMesh(tri);
        MeshGeodesy *mg = new MeshGeodesy(m_mesh);
        mg->compute_distances(source_vertex);
        mg->save_geodesic(output_file);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }

    return EXIT_SUCCESS;
}