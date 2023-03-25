#ifndef MESHGUI_INCLUDED
#define MESHGUI_INCLUDED

#include "MxGUI.hpp"

class MeshGUI : public MxGUI
{
public:
	MeshGUI();
	~MeshGUI();
	
	void initialize(int argc, char* argv[]);
	int add_menu_item(const char* name, int key, Fl_Callback *f, void* val=0, int flags=0);
};

#endif