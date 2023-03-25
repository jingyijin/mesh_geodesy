#include "MxGUI.hpp"
#include "MeshGUI.hpp"

MeshGUI::MeshGUI() 
{
}

MeshGUI::~MeshGUI() 
{}

void MeshGUI::initialize(int argc, char* argv[])
{
    MxGUI::initialize(argc, argv);
}

int MeshGUI::add_menu_item(const char* name, int key, Fl_Callback *f, void* val, int flags) 
{
    return MxGUI::menu_bar->add(name, key, f, val);
}
