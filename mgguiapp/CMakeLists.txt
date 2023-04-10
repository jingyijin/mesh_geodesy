cmake_minimum_required(VERSION 3.0)

project(mgguiapp)


find_package(FLTK REQUIRED)
if (NOT FLTK_FOUND)
    message(FATAL_ERROR "FLTK not found. Please install it using 'sudo apt-get install libfltk1.3-dev'.")
endif()

find_package(OpenGL REQUIRED)
if (NOT OpenGL_FOUND)
    message(FATAL_ERROR "OpenGL not found. Please install it using 'sudo apt-get install libgl1-mesa-dev'.")
endif()

find_library(XFT_LIBRARY NAMES Xft)
if (XFT_LIBRARY STREQUAL "XFT_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Xft library not found. Please install it using 'sudo apt-get install libxft-dev'.")
endif()

find_library(FONTCONFIG_LIBRARY NAMES fontconfig)
if (FONTCONFIG_LIBRARY STREQUAL "FONTCONFIG_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Fontconfig library not found. Please install it using 'sudo apt-get install libfontconfig1-dev'.")
endif()

find_library(XRENDER_LIBRARY NAMES Xrender)
if (XRENDER_LIBRARY STREQUAL "XRENDER_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Xrender library not found. Please install it using 'sudo apt-get install libxrender-dev'.")
endif()

find_library(XFIXES_LIBRARY NAMES Xfixes)
if (XFIXES_LIBRARY STREQUAL "XFIXES_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Xfixes library not found. Please install it using 'sudo apt-get install libxfixes-dev'.")
endif()

find_library(XCURSOR_LIBRARY NAMES Xcursor)
if (XCURSOR_LIBRARY STREQUAL "XCURSOR_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Xcursor library not found. Please install it using 'sudo apt-get install libxcursor-dev'.")
endif()

find_library(XINERAMA_LIBRARY NAMES Xinerama)
if (XINERAMA_LIBRARY STREQUAL "XINERAMA_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Xinerama library not found. Please install it using 'sudo apt-get install libxinerama-dev'.")
endif()

find_package(glog REQUIRED)
if (NOT glog_FOUND)
    message(FATAL_ERROR "glog not found. Please install it using 'sudo apt-get install libgoogle-glog-dev'.")
endif()


set(MGGUIAPP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(MGGUIAPP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(MGGUIAPP_SOURCES
    ${MGGUIAPP_SOURCE_DIR}/mxGUI.cpp
    ${MGGUIAPP_SOURCE_DIR}/meshGUI.cpp
    ${MGGUIAPP_SOURCE_DIR}/arcball.cpp
    ${MGGUIAPP_SOURCE_DIR}/baseball.cpp
    ${MGGUIAPP_SOURCE_DIR}/main.cpp
)

set(MGGUIAPP_HEADERS
    ${MGGUIAPP_INCLUDE_DIR}/mxGUI.hpp
    ${MGGUIAPP_INCLUDE_DIR}/meshGUI.hpp
    ${MGGUIAPP_INCLUDE_DIR}/arcball.hpp
    ${MGGUIAPP_INCLUDE_DIR}/baseball.hpp
)

set(CMAKE_BUILD_TYPE Debug)

add_executable(mgguiapp ${MGGUIAPP_SOURCES} ${MGGUIAPP_HEADERS})

target_include_directories(mgguiapp PRIVATE ${MGGUIAPP_INCLUDE_DIR})

target_link_libraries(mgguiapp mglib
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
