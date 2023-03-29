# version 3.11 or later of CMake or needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(mglib)

set(MGLIB_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(MGLIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(GFX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gfx/source")
set(GFX_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gfx/include")

set(MGLIB_SOURCES
    ${MGLIB_SOURCE_DIR}/general_math.cpp
    ${MGLIB_SOURCE_DIR}/geomesh.cpp
    ${MGLIB_SOURCE_DIR}/trimesh.cpp
    ${GFX_SOURCE_DIR}/heap.cpp
    ${GFX_SOURCE_DIR}/quat.cpp
)

set(MGLIB_HEADERS
    ${MGLIB_INCLUDE_DIR}/edgestruct.hpp
    ${MGLIB_INCLUDE_DIR}/general_math.hpp
    ${MGLIB_INCLUDE_DIR}/geomesh.hpp
    ${MGLIB_INCLUDE_DIR}/halfedge.hpp
    ${MGLIB_INCLUDE_DIR}/interval.hpp
    ${MGLIB_INCLUDE_DIR}/last_step.hpp
    ${MGLIB_INCLUDE_DIR}/manifold.hpp
    ${MGLIB_INCLUDE_DIR}/propagation.hpp
    ${MGLIB_INCLUDE_DIR}/trimesh.hpp
    ${MGLIB_INCLUDE_DIR}/ray.hpp
    ${MGLIB_INCLUDE_DIR}/typedef.hpp
    ${GFX_INCLUDE_DIR}/gfx.hpp
    ${GFX_INCLUDE_DIR}/geom3d.hpp
    ${GFX_INCLUDE_DIR}/heap.hpp
    ${GFX_INCLUDE_DIR}/mat3.hpp
    ${GFX_INCLUDE_DIR}/mat4.hpp
    ${GFX_INCLUDE_DIR}/quat.hpp
    ${GFX_INCLUDE_DIR}/vec2.hpp
    ${GFX_INCLUDE_DIR}/vec3.hpp
    ${GFX_INCLUDE_DIR}/vec4.hpp
)

add_library(mglib STATIC 
    ${MGLIB_SOURCES} 
    ${MGLIB_HEADERS}
)

target_include_directories(mglib PUBLIC 
    ${MGLIB_INCLUDE_DIR}
    ${GFX_INCLUDE_DIR}
)