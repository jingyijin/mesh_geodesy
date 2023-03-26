#include <iostream>
#include "MeshGUI.hpp"

using namespace std;

MeshGUI gui;

int main(int argc, char *argv[])
{
    gui.initialize(argc, argv);
    gui.toplevel->label("Mesh Geodesy");
    gui.load_mesh("../model/hand.obj");
    
    return gui.run();
}
