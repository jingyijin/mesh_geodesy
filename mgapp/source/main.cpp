#include <iostream>
#include "meshGUI.hpp"

using namespace std;

MeshGUI gui;

int main(int argc, char *argv[])
{
    gui.initialize(argc, argv);
    gui.m_toplevel->label("Mesh Geodesy");
    gui.load_mesh("../model/bunny.obj");

    return gui.run();
}
