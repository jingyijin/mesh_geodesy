#ifndef GFXGUI_INCLUDED // -*- C++ -*-
#define GFXGUI_INCLUDED

/************************************************************************

  Minimalist GUI framework.

  This package implements a baseline GUI framework for use in
  GFX-based applications.  Only a very specific kind of interface is
  supported: one where the application window consists primarily of an
  OpenGL drawing canvas.

  $Id: gui.h 443 2005-06-14 00:53:40Z garland $

 ************************************************************************/

#include <GL/gl.h>
#include <GL/glu.h>

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Output.H>

#include <string>

using namespace std;

class MxGUI;

enum { IMG_PNM=0, IMG_PNG=1, IMG_TIFF=2, IMG_JPEG=3, IMG_LIMIT=4 };

class MxGLCanvas : public Fl_Gl_Window
{
private:
    int m_last_click[2];
    MxGUI *m_app;

public:
    /**
     * @brief Overrides the draw method of Fl_Gl_Window to perform custom OpenGL rendering.
     */
    virtual void draw();
    /**
     * @brief Overrides the handle method of Fl_Gl_Window to handle events.
     * @param event The event to handle.
     * @return int Returns 1 if the event was handled, 0 otherwise.
     */
    virtual int handle(int event);
    /**
     * @brief Overrides the resize method of Fl_Gl_Window to handle canvas resizing.
     * @param x X position of the canvas.
     * @param y Y position of the canvas.
     * @param w Width of the canvas.
     * @param h Height of the canvas.
     */
    virtual void resize(int x, int y, int w, int h);

public:
    /**
    * @brief Default constructor for MxGLCanvas.
    * @param x X position of the canvas.
    * @param y Y position of the canvas.
    * @param w Width of the canvas.
    * @param h Height of the canvas.
    * @param label Label for the canvas.
    */
    MxGLCanvas(int x, int y, int w, int h, const char *label=NULL);
    /**
     * @brief Attaches a parent MxGUI object to this canvas.
     * @param a Pointer to the parent MxGUI object.
     */    
    void attach_app(MxGUI *a);
};

/**
* @class MxGUI
* @brief Class for creating a graphical user interface using FLTK.
* This class provides a simple interface for creating a graphical user interface (GUI)
* using the Fast Light Toolkit (FLTK). It provides a toplevel window, a GL canvas for
* rendering 3D graphics, a status line, and a menu bar. The GUI can be initialized by
* calling the initialize method and then run by calling the run method. The animate
* method can be used to toggle animation on or off. The status method can be used to
* display a status message in the status line.
*/
class MxGUI
{
private:
    int m_w_offset, m_h_offset;
    Fl_Window *create_window(int xw=640, int yw=480, int pad=5);

public:
    /**
    * @brief Pointer to the top-level FLTK window.
    */
    Fl_Window *m_toplevel;
    /**
     * @brief Pointer to the OpenGL canvas for rendering 3D graphics.
     */    
    MxGLCanvas *m_canvas;
    /**
     * @brief Pointer to the status line.
     */    
    Fl_Output *m_status_line;
    /**
     * @brief Pointer to the menu bar.
     */    
    Fl_Menu_Bar *m_menu_bar;
    /**
     * @brief Pointer to the layout menu item.
     */    
    Fl_Menu_Item *m_menu_layout;
    /**
     * @brief Default frames per second for animation.
     */
    float m_default_fps;
    /**
     * @brief Target frames per second for animation.
     */    
    float m_target_fps;

    /**
     * @brief Pointer to the current MxGUI object.
     */
    static MxGUI *m_current;	// There should only be one.

    /**
     * @brief Constructor for MxGUI.
     */
    MxGUI();
    /**
     * @brief Virtual destructor for MxGUI.
     */    
    virtual ~MxGUI() {}

    /**
     * @brief Initializes the GUI.
     * @param argc The number of command line arguments.
     * @param argv The command line arguments.
     * @param layout Pointer to the layout menu item.
     * @param xw The width of the window.
     * @param yw The height of the window.
     */
    virtual void initialize(int argc, char **argv,
                Fl_Menu_Item *layout=NULL,
                int xw=640, int yw=480);
    /**
     * @brief Runs the GUI.
     * @return An integer representing the exit status.
     */                
    virtual int run();

    /**
     * @brief Displays a status message in the status line.
     * @param fmt The format string for the message.
     * @return An integer representing the number of characters printed.
     */
    int status(const char *fmt, ...);
    /**
     * @brief Toggles animation on or off.
     * @param will A boolean value indicating whether animation will be turned on or off.
     */
    void animate(bool will);
//    bool snapshot_to_file(int format, const char *filenamep=NULL);
    void resize_canvas(int width, int height);
    void lock_size();
    void unlock_size();

    void title(const char *l) { m_toplevel->label(l); }

    // Menu construction and standard callbacks
    int add_menu(const string&, int key, Fl_Callback *cb, int flags=0);
    int add_toggle_menu(const string&, int key, bool& val, int flags=0);
    static void cb_toggle(Fl_Menu_ *m, bool *flag);

public:
    //
    // Callback functions that get executed in response to menu commands.
    virtual void cb_new();
    virtual void cb_exit();
    virtual void cb_snapshot(int);
    virtual void cb_animate(Fl_Menu_ *m);
    virtual void cb_fps();
    virtual void cb_vga_size(int width);  // uses 4:3  aspect ratio
    virtual void cb_hdtv_size(int width); // uses 16:9 aspect ratio
    virtual void cb_dv_size(int width);   // uses 3:2  aspect ratio

    virtual void cb_save_view_to_file();
    virtual void cb_load_view_from_file();
    virtual bool save_view_to_file();
    virtual bool load_view_from_file();

public:
    //
    // Applications are customized by overriding the following methods.

    // Override these methods to control the contents of the GL canvas
    virtual void setup_for_drawing();
    virtual void draw_contents();
    virtual void update_animation();

    // Override these methods to receive events from the GL canvas
    virtual bool mouse_down(int *where, int which);
    virtual bool mouse_up(int *where, int which);
    virtual bool mouse_drag(int *where, int *last, int which);
    virtual bool key_press(int key);

    // Override these methods to get command line arguments
    virtual int cmdline_option(int argc, char **argv, int& index);
    virtual void cmdline_file(const char *file);

    // Override these methods to add custom interface elements
    virtual void add_upper_controls(int& yfill, const int pad) {}
    virtual void add_lower_controls(int& yfill, const int pad) {}

    // Override this method to free memory, close files, etc.
    virtual void cleanup_for_exit() {}
};

////////////////////////////////////////////////////////////////////////
//
// This template makes it easier to create FLTK-compliant callbacks.
// In particular, its purpose is to construct static thunks for
// calling member functions of MxGUI-derived classes.
//

template<class Gui>
struct MxBinder
{
    typedef void (Gui::*GuiCommand)();
    typedef void (Gui::*GuiCommand1)(int);
    typedef void (Gui::*GuiCommand2)(Fl_Menu_ *);

    template<GuiCommand cmd>
    static void to(Fl_Widget *, void *data)
    {
    Gui *gui = static_cast<Gui*>(data);
    (gui->*cmd)();
    gui->m_canvas->redraw();
    }

    template<GuiCommand2 cmd>
    static void to_menu(Fl_Widget *w, void *data)
    {
    Gui *gui = static_cast<Gui*>(data);
    (gui->*cmd)(static_cast<Fl_Menu_ *>(w));
    gui->m_canvas->redraw();
    }

    template<GuiCommand1 cmd, int i>
    static void to_arg(Fl_Widget *, void *data)
    {
    Gui *gui = static_cast<Gui*>(data);
    (gui->*cmd)(i);
    gui->m_canvas->redraw();
    }
};

////////////////////////////////////////////////////////////////////////
//
// These macros make static FLTK menu definitions look a little nicer.
//

#define MXGUI_BEGIN_MENU(name) {name, 0, 0, 0, FL_SUBMENU},
#define MXGUI_END_MENU {0},
#define MXGUI_FINISH_MENUBAR {0}

// GFXGUI_INCLUDED
#endif
