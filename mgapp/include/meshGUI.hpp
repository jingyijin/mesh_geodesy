#ifndef MESHGUI_INCLUDED
#define MESHGUI_INCLUDED

#include "mxGUI.hpp"
#include "trimesh.hpp"
#include "arcball.hpp"
#include "vec4.hpp"

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
        glMaterialfv(FACE, GL_EMISSION, emit);
        glMaterialfv(FACE, GL_AMBIENT, r_amb);
        glMaterialf(FACE, GL_SHININESS, shininess);
        new_diff = rgbRefl(0.9f, 0.25f, 0.3f, 0.5f);;
        glMaterialfv(FACE, GL_DIFFUSE, new_diff);
        new_diff = 0.3f * new_diff;
        glMaterialfv(FACE, GL_SPECULAR, r_spec);
    }
};

class MeshGUI : public MxGUI
{
protected:
    Arcball ball;
    Vec3 bb_min, bb_max;

    Material* mat;

public:
    TriMesh *mesh;
    enum {draw_mode_wireframe, draw_mode_color, draw_mode_solid} draw_mode;
    enum {NoSelect, FSelect, VSelect} selection_mode;

    // Drawing protocol
    bool will_draw_bbox;
    bool will_draw_surface_fnormal;
    bool will_draw_mesh;

public:
    MeshGUI();
    ~MeshGUI();
    
    void initialize(int argc, char* argv[]);
    int add_menu_item(const char* name, int key, Fl_Callback *f, void* val=0, int flags=0);

    void setup_for_drawing(); 
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

    bool mouse_down(int *where, int which); 
    bool mouse_drag(int *where, int *last, int which); 
    bool mouse_up(int *where, int which);
    bool key_press(int key); 

    // callback functions
    static void cb_open_file();
    static void cb_save_file();
    void load_mesh(const string& filename);
};

extern MeshGUI gui;

#endif