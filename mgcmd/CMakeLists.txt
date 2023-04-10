cmake_minimum_required(VERSION 3.0)

project(mgcmd)

find_package(glog REQUIRED)
if (NOT glog_FOUND)
    message(FATAL_ERROR "glog not found. Please install it using 'sudo apt-get install libgoogle-glog-dev'.")
endif()


set(MGCMD_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/source")
set(MGCMD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(MGCMD_SOURCES
    ${MGCMD_SOURCE_DIR}/main.cpp
)

set(MGCMD_HEADERS
)

set(CMAKE_BUILD_TYPE Debug)

add_executable(mgcmd ${MGCMD_SOURCES} ${MGCMD_HEADERS})

target_include_directories(mgcmd PRIVATE ${MGCMD_INCLUDE_DIR})

target_link_libraries(mgcmd mglib
    glog::glog
)
