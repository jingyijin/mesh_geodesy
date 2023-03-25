#include <iostream>
#include "MeshGUI.hpp"

using namespace std;

MeshGUI gui;

int main(int argc, char *argv[])
{
	gui.initialize(argc, argv);
	gui.resize_canvas(800, 500);
    gui.title("Mesh Geodesy");
	cout << "signs of life" << endl;
    return gui.run();
}
