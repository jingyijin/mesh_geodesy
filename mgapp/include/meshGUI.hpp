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

/**
* @brief A class for materials with Phong illumination model.
* This class defines the properties of a material used in Phong illumination model.
* It defines the ambient, diffuse, and specular reflection properties, the amount of
* light emitted by the material, the amount of light reflected by the material, and
* the shininess exponent of the material.
*/
class Material
{
public:
    rgbColor emit;      // Light emitted by this material

    rgbRefl r_amb;      // Fraction of ambient light reflected
    rgbRefl r_diff;     // Fraction of incoming light reflected diffusely 
    rgbRefl r_spec;     // Fraction of incoming light reflected specularly
    rgbRefl r_transp_diff;     // Fraction of incoming light reflected diffusely 
    double shininess;   // Exponent for Phong illumination model

    /**
    * @brief Construct a new Material object with default properties.
    * The default material properties define a bronze-like material.
    */
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

    /**
    * @brief Load standard material properties into the OpenGL pipeline.
    * This method loads the ambient, diffuse, and specular reflection properties,
    * and the shininess exponent of the material into the OpenGL pipeline for
    * both front and back faces.
    */
    void load_standard() const 
    {
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, r_amb);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, r_diff);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, r_spec);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
    }

    /**
    * @brief Load yellow material properties into the OpenGL pipeline.
    */
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

    /**
    * @brief Load red material properties into the OpenGL pipeline.
    */
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

/**
* @brief The MeshGUI class is a subclass of MxGUI that provides an interface for visualizing and interacting with 3D meshes.
*/
class MeshGUI : public MxGUI
{
protected:
    Arcball m_ball;            // Arcball for mouse interaction
    Vec3 m_bb_min, m_bb_max;   // Bounding box of the mesh

    Material *m_mat;           // Material for rendering the mesh
    GLUquadricObj *m_obj;      // Quadric object for rendering the selected vertex

    int m_grid_period;         // Period of the texture grid
public:
    GeoTriMesh *m_mesh;        // The core mesh data structure
    MLS *m_mls;                // The data structure to trigger distance calculation
    ByteRaster* m_texture;     // The texture image

    enum {Draw_mode_wireframe, Draw_mode_color, Draw_mode_solid} m_draw_mode;   // Drawing mode
    enum {Noselect, Fselect, Vselect} m_selection_mode;                         // Selection mode

    // Drawing protocol
    bool m_will_draw_bbox;                  // Flag to draw bounding box
    bool m_will_draw_surface_fnormal;       // Flag to draw surface face normals
    bool m_will_draw_mesh;                  // Flag to draw mesh in wireframe
    bool m_will_draw_vertices;              // Flag to draw vertices
    bool m_will_draw_geodesic_distance;     // Flag to draw geodesic distance
    bool m_will_draw_geodesic_path;         // Flag to draw geodesic path

    int m_selected_vertex;                  // Selected vertex

    // Animation
    int m_displayed_face_id;            // The cutoff id for face to be displayed in the animation
    int m_display_face_increment;       // The increment of the face ids to be displayed in the next frame of the animation

public:
    /**
    * @brief MeshGUI constructor.
    */
    MeshGUI();
    /**
     * @brief MeshGUI destructor.
     */    
    ~MeshGUI();
    
    /**
     * @brief Initializes the MeshGUI object.
     *
     * @param argc Number of command line arguments.
     * @param argv Array of command line arguments.
     */
    void initialize(int argc, char* argv[]);
    /**
     * @brief Adds a menu item to the GUI.
     *
     * @param name Name of the menu item.
     * @param key Key associated with the menu item.
     * @param f Callback function to be executed when the menu item is selected.
     * @param val Pointer to the value associated with the menu item.
     * @param flags Flags associated with the menu item.
     *
     * @return The ID of the newly added menu item.
     */    
    int add_menu_item(const char* name, int key, Fl_Callback *f, void* val=0, int flags=0);

    /**
     * @brief Sets up the MeshGUI object for drawing.
     */
    void setup_for_drawing(); 
    /**
     * @brief Sets up the texture for the MeshGUI object.
     */
    void setup_texture();
    /**
     * @brief Sets up the default texture for the MeshGUI object.
     */    
    void default_texture();
    /**
     * @brief Begins redrawing the MeshGUI object.
     */
    void begin_redraw();
    /**
     * @brief Sets up the default redrawing for the MeshGUI object.
     */
    void default_redraw();
    /**
     * @brief Ends redraw for the MeshGUI object.
     */
    void end_redraw();

    /**
    * @brief Apply the current camera transformation to the OpenGL context.
    */
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
    void update_animation();
    void animate(bool will);

    int pick_vertex(int where[2]); 
    void up_frequency();
    void down_frequency();

    // callback functions
    static void cb_open_file();
    static void cb_save_file();    void apply_camera();
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
    void update_animation();
    void animate(bool will);

    int pick_vertex(int where[2]); 
    void up_frequency();
    void down_frequency();

    // callback functions
    static void cb_open_file();
    static void cb_save_file();
    static void cb_load_default_texture();
    static void cb_load_texture();
    static void cb_save_distance();
    static void cb_load_distance();

    void load_mesh(const string& filename);
};
extern MeshGUI gui;

#endif