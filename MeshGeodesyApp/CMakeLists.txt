cmake_minimum_required(VERSION 3.0)

project(MeshGeodesyApp)

find_package(FLTK REQUIRED)
find_package(OpenGL REQUIRED)
find_library(XFT_LIBRARY NAMES Xft)
find_library(FONTCONFIG_LIBRARY NAMES fontconfig)
find_library(XRENDER_LIBRARY NAMES Xrender)
find_library(XFIXES_LIBRARY NAMES Xfixes)
find_library(XCURSOR_LIBRARY NAMES Xcursor)
find_library(XINERAMA_LIBRARY NAMES Xinerama)

set(MESHAPP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(MESHAPP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(MESHAPP_SOURCES
    ${MESHAPP_SOURCE_DIR}/MxGUI.cpp
    ${MESHAPP_SOURCE_DIR}/MeshGUI.cpp
    ${MESHAPP_SOURCE_DIR}/main.cpp
)

set(MESHAPP_HEADERS
    ${MESHAPP_INCLUDE_DIR}/MxGUI.hpp
    ${MESHAPP_INCLUDE_DIR}/MeshGUI.hpp
)

add_executable(MeshGeodesyApp ${MESHAPP_SOURCES} ${MESHAPP_HEADERS})

target_include_directories(MeshGeodesyApp PRIVATE ${MESHAPP_INCLUDE_DIR})

target_link_libraries(MeshGeodesyApp MeshGeodesyLib
    ${FLTK_LIBRARIES}
    ${OPENGL_LIBRARIES} 
    ${XFT_LIBRARY}
    ${FONTCONFIG_LIBRARY}
    ${XRENDER_LIBRARY}
    ${XFIXES_LIBRARY}
    ${XCURSOR_LIBRARY}
    ${XINERAMA_LIBRARY}
)
