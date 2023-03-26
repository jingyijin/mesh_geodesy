
#include "MxGUI.hpp"
#include "MeshGUI.hpp"
#include "vec4.hpp"

#include <FL/Fl_File_Chooser.H>

MeshGUI::MeshGUI() : MxGUI(), 
		draw_mode(draw_mode_wireframe), 
		will_draw_bbox(false), 
		will_draw_surface_fnormal(true), 
		will_draw_mesh(true)	
{
	mesh = new TriMesh();
	mat = new Material();
}

MeshGUI::~MeshGUI() 
{
	if (mesh) delete mesh;
	if (mat)  delete mat;
}

void MeshGUI::initialize(int argc, char* argv[])
{
    MxGUI::initialize(argc, argv);

	cout << "MeshGUI::initialize()" << endl;
	
	// Set up the menu bar
	add_menu_item("&File/Open", FL_CTRL+'o', (Fl_Callback *)cb_open_file, NULL);
	add_menu_item("&File/Open", FL_CTRL+'o', (Fl_Callback *)cb_open_file, NULL);
	add_menu_item("&File/Save", FL_CTRL+'s', (Fl_Callback *)cb_save_file, NULL);
}

int MeshGUI::add_menu_item(const char* name, int key, Fl_Callback *f, void* val, int flags) 
{
    return MxGUI::menu_bar->add(name, key, f, val);
}

void MeshGUI::cb_open_file() 
{
	string input_filename = fl_file_chooser("Select input file", "*{.obj}", NULL);

	if (!input_filename.empty())
		gui.load_mesh(input_filename);
}

void MeshGUI::cb_save_file()
{
	string output_filename = fl_file_chooser("Select output file", "*{.obj}", NULL);

 	if (!output_filename.empty()) {
		if (gui.mesh)
			gui.mesh->write_to_file(output_filename);
	}
}

void MeshGUI::load_mesh(const string& filename)
{
	if (mesh) mesh->read_from_file(filename);
	
	mesh->initialize();
	reset_camera();
	canvas->redraw();
}

void MeshGUI::setup_for_drawing() 
{
	cout << "MeshGUI::setup_for_drawing()" << endl;

    glClearColor(1.f, 1.f, 1.f, 0.0f);
	
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
	glEnableClientState(GL_VERTEX_ARRAY);
	
    // Enable lighting and set up the lighting environment.  We specify a
    // global ambient glow and create two point lights.
    glEnable(GL_LIGHTING);
    Vec4f ambient_light(1.0, 1.0, 1.0, 1.0);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, (float *)ambient_light);
	
    const Vec4f light0_pos(0.0f, 0.5f, 1.0f, 0.0f);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glEnable(GL_LIGHT0);

	mat->load_standard();

    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

}

void MeshGUI::draw_contents()
{
	cout << "MeshGUI::draw_contents()" << endl;
    begin_redraw();

	glPushMatrix();
    default_redraw();
	glPopMatrix();
    
	end_redraw();
}

void MeshGUI::begin_redraw() 
{
	cout << "MeshGUI::begin_redraw()" << endl;
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
	if( selection_mode == NoSelect ) 
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
	}
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	apply_camera();
	
	ball.apply_transform();
}

void MeshGUI::default_redraw() 
{
	if (!mesh) return;
	
	draw_surface_fnormal();
}

void MeshGUI::end_redraw() 
{
	ball.unapply_transform();
}

void MeshGUI::apply_camera() 
{
	float aspect = (float)canvas->w() / (float)canvas->h();
	camera_lookat(bb_min, bb_max, aspect);
}

void MeshGUI::reset_camera() 
{
	if (mesh) mesh->computeBBox(bb_min, bb_max);
	Vec3 ctr = (bb_max + bb_min)/2.0;
	double radius = norm(bb_max - ctr);
	ball.bounding_sphere(ctr, radius);
}

void MeshGUI::camera_lookat(const Vec3& min, const Vec3& max, double aspect)
{
    Vec3 up(0, 1, 0);
    double fovy = 60.0;

    Vec3 at = (max + min)/2.0;         // look at the center of the bbox
    double radius = norm(max - at);    // radius of a bounding sphere
    double d = 3*radius / tan(fovy * M_PI/180.0);

    Vec3 from = at;
    from[2] += d;

    double znear = d/20;
    double zfar = 10*d;

    glMatrixMode(GL_PROJECTION);
    gluPerspective(fovy, aspect, znear, zfar);


    glMatrixMode(GL_MODELVIEW);
    gluLookAt(from[0], from[1], from[2],
	      at[0], at[1], at[2],
	      up[0], up[1], up[2]);
}

void MeshGUI::draw_mesh() 
{
	glDisable(GL_POLYGON_OFFSET_FILL);
	
	glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
	glVertexPointer(3, GL_DOUBLE, 0, &mesh->vertex[0]);
	
	glDisable(GL_LIGHTING);
	glColor3f(0.f, 0.f, 0.f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	TriMesh::FaceList& face = mesh->face;
	int fsize = face.size();

	glBegin(GL_TRIANGLES);
	for (int i=0; i<fsize; i++) 
	{
		const Face& f = face[i];
		
		setup_face_state(i);
		glArrayElement(f[0]);
		glArrayElement(f[1]);
		glArrayElement(f[2]);
	}
	glEnd();
	
	if (!will_draw_surface_fnormal) 
	{
		GLfloat bkg_color[4];
		glGetFloatv(GL_COLOR_CLEAR_VALUE, bkg_color);
		glColor4fv(bkg_color);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		draw_surface_fnormal();
	}
	
	glPopAttrib();
}

void MeshGUI::draw_surface_fnormal() 
{
	glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
	glVertexPointer(3, GL_DOUBLE, 0, &mesh->vertex[0]);
	
	if (will_draw_mesh) 
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}
	
	if (mesh->fnormal.empty())
		mesh->computeFNormal();

	TriMesh::FaceList& face = mesh->face;
	int fsize = face.size();
	if (draw_mode == draw_mode_solid) 
	{
		glBegin(GL_TRIANGLES);
		for (int i=0; i<fsize; i++) 
		{
			const Face& f = face[i];
			
			setup_face_state(i); 
			glArrayElement(f[0]);
			glArrayElement(f[1]);
			glArrayElement(f[2]);
		}
		glEnd();
	} 
	else if (draw_mode == draw_mode_wireframe) 
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glColor3f(0.f, 0.f, 0.f);
		glDisable(GL_LIGHTING);
		glBegin(GL_TRIANGLES);
		for (int i=0; i<fsize; i++) 
		{
			const Face& f = face[i];
			
			setup_face_state(i); 
			glArrayElement(f[0]);
			glArrayElement(f[1]);
			glArrayElement(f[2]);
		}
		glEnd();
	} 
	glPopAttrib();
}

void MeshGUI::setup_face_state(int fid) 
{
    if (fid < mesh->fnormal.size())
		glNormal3dv(mesh->fnormal[fid]);
}

bool MeshGUI::mouse_down(int *where, int which) 
{
	cout << "mouse_down" << endl;
    return ball.mouse_down(where, which);
}

bool MeshGUI::mouse_drag(int *where, int *last, int which) 
{
	cout << "mouse_drag" << endl;
    if (which == 1 && Fl::event_state(FL_ALT))
		return ball.mouse_drag(where, last, 3);
    else if (which == 1 && Fl::event_state(FL_SHIFT))
		return ball.mouse_drag(where, last, 2);
    return ball.mouse_drag(where, last, which);
}

bool MeshGUI::mouse_up(int *where, int which) 
{
	cout << "mouse_up" << endl;
    return ball.mouse_up(where, which);
}