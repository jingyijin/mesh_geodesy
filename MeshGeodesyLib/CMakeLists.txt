# version 3.11 or later of CMake or needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(MeshGeodesyLib)

set(MESHLIB_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(MESHLIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(MESHLIB_SOURCES
    ${MESHLIB_SOURCE_DIR}/trimesh.cpp
)

set(MESHLIB_HEADERS
    ${MESHLIB_INCLUDE_DIR}/trimesh.hpp
)

add_library(MeshGeodesyLib STATIC ${MESHLIB_SOURCES} ${MESHLIB_HEADERS})

target_include_directories(MeshGeodesyLib PUBLIC ${MESHLIB_INCLUDE_DIR})