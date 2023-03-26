# version 3.11 or later of CMake or needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(MeshGeodesyLib)

set(MESHGEODESYLIB_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(MESHGEODESYLIB_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(GFX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gfx/source")
set(GFX_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gfx/include")

set(MESHGEODESYAPP_SOURCES
    ${MESHGEODESYLIB_SOURCE_DIR}/trimesh.cpp
    ${GFX_SOURCE_DIR}/quat.cpp
)

set(MESHGEODESYAPP_HEADERS
    ${MESHGEODESYLIB_INCLUDE_DIR}/trimesh.hpp
    ${GFX_INCLUDE_DIR}/gfx.hpp
    ${GFX_INCLUDE_DIR}/geom3d.hpp
    ${GFX_INCLUDE_DIR}/mat3.hpp
    ${GFX_INCLUDE_DIR}/mat4.hpp
    ${GFX_INCLUDE_DIR}/quat.hpp
    ${GFX_INCLUDE_DIR}/vec2.hpp
    ${GFX_INCLUDE_DIR}/vec3.hpp
    ${GFX_INCLUDE_DIR}/vec4.hpp
)

add_library(MeshGeodesyLib STATIC 
    ${MESHGEODESYAPP_SOURCES} 
    ${MESHGEODESYAPP_HEADERS}
)

target_include_directories(MeshGeodesyLib PUBLIC 
    ${MESHGEODESYLIB_INCLUDE_DIR}
    ${GFX_INCLUDE_DIR}
)