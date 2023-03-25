# version 3.11 or later of CMake needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(
    MeshGeodesy
    VERSION 0.1
    DESCRIPTION "Mesh Geodesy - an intuitive program that computes \
        geodesic distance on triangle meshes.")

add_subdirectory(MeshGeodesyApp)      # look in app subdirectory for CMakeLists.txt to process
add_subdirectory(MeshGeodesyLib)      # look in lib subdirectory for CMakeLists.txt to process
