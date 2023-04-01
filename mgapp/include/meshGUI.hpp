#ifndef MESHGUI_INCLUDED
#define MESHGUI_INCLUDED

/************************************************************************
 * MeshGUI is a class that provides a GUI for a mesh.  It is a subclass
 * of mxGUI, which provides the basic GUI functionality.  The class
 * provides a number of virtual functions that can be overridden to
 * provide custom functionality.  The class provides customized rendering
 * of the mesh, and provides a number of callbacks for mouse and keyboard
 * events.
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "mxGUI.hpp"
#include "geomesh.hpp"
#include "arcball.hpp"
#include "vec4.hpp"
#include "mls.hpp"
#include "raster.hpp"

typedef Vec4f rgbColor;
typedef Vec4f rgbRefl;

class Material
{
public:
    rgbColor emit;      // Light emitted by this material

    rgbRefl r_amb;      // Fraction of ambient light reflected
    rgbRefl r_diff;     // Fraction of incoming light reflected diffusely 
    rgbRefl r_spec;     // Fraction of incoming light reflected specularly
    rgbRefl r_transp_diff;     // Fraction of incoming light reflected diffusely 
    double shininess;   // Exponent for Phong illumination model

    Material() 
    {
        Vec4f rgb = Vec4f(0.912f, 0.717f, 0.505f, 1.0f);
        Vec4f rgb_transp = Vec4f(0.9f, 0.9f, 0.3f, 0.3f);
        r_amb  = 0.1f*rgb;
        r_diff = 1.0f*rgb;
        r_spec = 0.3f*rgb;
        r_transp_diff = rgb_transp;
        shininess = 100;
    }

    void load_standard() const 
    {
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, r_amb);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, r_diff);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, r_spec);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
    }

    void load_yellow() const 
    {
        const GLenum FACE = GL_FRONT;

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

        glMaterialfv(FACE, GL_EMISSION, emit);
        glMaterialfv(FACE, GL_AMBIENT, r_amb);
        glMaterialf(FACE, GL_SHININESS, shininess);
        glMaterialfv(FACE, GL_DIFFUSE, r_transp_diff);
        glMaterialfv(FACE, GL_SPECULAR, r_spec);

    }

    void load_red() const
    {
        const GLenum FACE = GL_FRONT;

        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

        rgbRefl new_diff;
        
        // Copy this material's values into the current OpenGL material
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emit);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, r_amb);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
        new_diff = rgbRefl(0.9f, 0.1f, 0.1f, 1.0f);;
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, new_diff);
        new_diff = 0.3f * new_diff;
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, r_spec);
    }
};

class MeshGUI : public MxGUI
{
protected:
    Arcball m_ball;
    Vec3 m_bb_min, m_bb_max;

    Material *m_mat;
    GLUquadricObj *m_obj;

    int m_grid_period;
public:
    GeoTriMesh *m_mesh;
    MLS *m_mls;
   	ByteRaster* m_texture;

    enum {Draw_mode_wireframe, Draw_mode_color, Draw_mode_solid} m_draw_mode;
    enum {Noselect, Fselect, Vselect} m_selection_mode;

    // Drawing protocol
    bool m_will_draw_bbox;
    bool m_will_draw_surface_fnormal;
    bool m_will_draw_mesh;
    bool m_will_draw_vertices;
	bool m_will_draw_geodesic_distance;
	bool m_will_draw_geodesic_path;

    int m_selected_vertex;

public:
    MeshGUI();
    ~MeshGUI();
    
    void initialize(int argc, char* argv[]);
    int add_menu_item(const char* name, int key, Fl_Callback *f, void* val=0, int flags=0);

    void setup_for_drawing(); 
    void setup_texture();
    void default_texture();
    void begin_redraw();
    void default_redraw();
    void end_redraw();

    void apply_camera();
    void reset_camera();
    void camera_lookat(const Vec3& min, const Vec3& max, double aspect);

    void setup_face_state(int fid);
    void draw_mesh();
    void draw_contents();
    void draw_bbox();
    void draw_box(const Vec3f& min, const Vec3f& max);
    void draw_surface_fnormal();
    void draw_for_selection();
    void draw_selection();
    void draw_vertices(); 
    void draw_geodesic_distance();
    void draw_geodesic_path();

    bool mouse_down(int *where, int which); 
    bool mouse_drag(int *where, int *last, int which); 
    bool mouse_up(int *where, int which);
    bool key_press(int key); 

    int pick_vertex(int where[2]); 
    void up_frequency();
    void down_frequency();

    // callback functions
    static void cb_open_file();
    static void cb_save_file();
    void load_mesh(const string& filename);
};

extern MeshGUI gui;

#endif