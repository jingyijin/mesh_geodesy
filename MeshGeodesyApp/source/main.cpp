#include <iostream>
#include "MeshGUI.hpp"

using namespace std;

MeshGUI gui;

int main(int argc, char *argv[])
{
	gui.initialize(argc, argv);
	gui.resize_canvas(800, 500);
    gui.title("Mesh Geodesy");
	gui.mesh->read_from_file("../model/hand.obj");
	gui.canvas->redraw();
	
    return gui.run();
}
