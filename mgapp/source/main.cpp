#include <iostream>
#include <glog/logging.h>

#include "meshGUI.hpp"

using namespace std;

MeshGUI gui;

int main(int argc, char *argv[])
{
    google::InitGoogleLogging(argv[0]);

    gui.initialize(argc, argv);
    gui.m_toplevel->label("Mesh Geodesy");
    gui.load_mesh("../model/model.obj");

    return gui.run();
}
