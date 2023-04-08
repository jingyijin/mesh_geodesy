cmake_minimum_required(VERSION 3.0)

project(mgapp)

find_package(FLTK REQUIRED)
find_package(OpenGL REQUIRED)
find_library(XFT_LIBRARY NAMES Xft)
find_library(FONTCONFIG_LIBRARY NAMES fontconfig)
find_library(XRENDER_LIBRARY NAMES Xrender)
find_library(XFIXES_LIBRARY NAMES Xfixes)
find_library(XCURSOR_LIBRARY NAMES Xcursor)
find_library(XINERAMA_LIBRARY NAMES Xinerama)
find_package(glog REQUIRED)

set(MGAPP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(MGAPP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(MGAPP_SOURCES
    ${MGAPP_SOURCE_DIR}/mxGUI.cpp
    ${MGAPP_SOURCE_DIR}/meshGUI.cpp
    ${MGAPP_SOURCE_DIR}/arcball.cpp
    ${MGAPP_SOURCE_DIR}/baseball.cpp
    ${MGAPP_SOURCE_DIR}/main.cpp
)

set(MGAPP_HEADERS
    ${MGAPP_INCLUDE_DIR}/mxGUI.hpp
    ${MGAPP_INCLUDE_DIR}/meshGUI.hpp
    ${MGAPP_INCLUDE_DIR}/arcball.hpp
    ${MGAPP_INCLUDE_DIR}/baseball.hpp
)

set(CMAKE_BUILD_TYPE Debug)

add_executable(mgapp ${MGAPP_SOURCES} ${MGAPP_HEADERS})

target_include_directories(mgapp PRIVATE ${MGAPP_INCLUDE_DIR})

target_link_libraries(mgapp mglib
    ${FLTK_LIBRARIES}
    ${OPENGL_LIBRARIES} 
    ${XFT_LIBRARY}
    ${FONTCONFIG_LIBRARY}
    ${XRENDER_LIBRARY}
    ${XFIXES_LIBRARY}
    ${XCURSOR_LIBRARY}
    ${XINERAMA_LIBRARY}
    glog::glog
)
