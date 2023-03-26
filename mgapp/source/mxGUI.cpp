/************************************************************************

  MxGUI
  
  $Id: gui.cxx 447 2005-06-18 14:01:16Z garland $

 ************************************************************************/

// #include <gfx/gfx.h>

#define HAVE_FLTK
#ifdef HAVE_FLTK

#include "mxGUI.hpp"

#include <FL/Fl_File_Chooser.H>
#include <FL/fl_ask.H>

#include <cstdio>
#include <cstdarg>
#include <fstream>
#include <iostream>

MxGUI *MxGUI::m_current = NULL;

MxGLCanvas::MxGLCanvas(int x, int y, int w, int h, const char *l)
    : Fl_Gl_Window(x, y, w, h, l)
{
    m_last_click[0] = m_last_click[1] = -1;
    m_app = NULL;
}

void MxGLCanvas::attach_app(MxGUI *a)
{
    if( !m_app )
    m_app = a;
}

void MxGLCanvas::resize(int x, int y, int w, int h)
{
    Fl_Gl_Window::resize(x, y, w, h);

    if( shown() )
    {
        make_current();
        glViewport(0, 0, w, h);
        invalidate();
    }
}

void MxGLCanvas::draw()
{
    if( !valid() )
    {
        valid(1);
        if(m_app) m_app->setup_for_drawing();
    }

    if(m_app) m_app->draw_contents();
}

int MxGLCanvas::handle(int event)
{
    bool need_redraw = false;

    int where[2];  where[0] = Fl::event_x();  where[1] = Fl::event_y();

    // NOTE: Normally, we examine event_state() rather than
    // event_button() because it is valid no matter what the generating
    // event whereas event_button() is only valid for PUSH and RELEASE
    // events.  However, since event_state() only tells us what buttons
    // are *pressed*, we need to revert to event_button() on RELEASE
    // events.
    //
    int which = Fl::event_button();

    if( event != FL_RELEASE )
    {
    if( Fl::event_state(FL_BUTTON1) )
    {
        // emulate middle button by combination of left & right buttons
        if( Fl::event_state(FL_BUTTON3) )  which = 2;
        else                               which = 1;
    }
    else if( Fl::event_state(FL_BUTTON2) ) which = 2;
    else if( Fl::event_state(FL_BUTTON3) ) which = 3;
    }

    switch( event )
    {
    case FL_PUSH:
    need_redraw = m_app && m_app->mouse_down(where, which);
    m_last_click[0]=where[0];  m_last_click[1]=where[1];
    break;

    case FL_RELEASE:
    need_redraw = m_app && m_app->mouse_up(where, which);
    break;

    case FL_DRAG:
    need_redraw = m_app && m_app->mouse_drag(where, m_last_click, which);
    m_last_click[0]=where[0];  m_last_click[1]=where[1];
    break;

    case FL_FOCUS:
    case FL_UNFOCUS:
    // Do nothing special
    break;

    case FL_KEYBOARD:
    if( !m_app || !m_app->key_press(Fl::event_key()) )
        return 0;
    break;

    default:
    return Fl_Gl_Window::handle(event);
    }

    if( need_redraw )
    redraw();

    return 1;
}

////////////////////////////////////////////////////////////////////////
//
// Default menu system.  Although the application can completely
// override this, it's expected that most programs will just add their
// own additional entries.
//

void MxGUI::cb_new() { }

void MxGUI::cb_exit() { cleanup_for_exit(); exit(0); }

void MxGUI::cb_animate(Fl_Menu_ *m)
    { animate(m->mvalue()->value()!=0); }

void MxGUI::cb_fps()
{
    // Convert default_fps to a string
    static char fps[64]; sprintf(fps, "%.1f", m_default_fps);

    const char *result = fl_input("Number of frames per second to draw", fps);
    if( result )
    {
    m_default_fps = atof(result);
    if( m_target_fps>0 ) m_target_fps=m_default_fps;
    }
}

void MxGUI::cb_snapshot(int format)
{
    m_canvas->redraw();		// don't want to snap menu
//    snapshot_to_file(format);	// snapshot what's drawn
}

void MxGUI::cb_vga_size(int xw)
{
    if( m_toplevel->resizable() )
    resize_canvas(xw, (3*xw)/4);
}

void MxGUI::cb_hdtv_size(int xw)
{
    if( m_toplevel->resizable() )
    MxGUI::m_current->resize_canvas(xw, (9*xw)/16);
}

void MxGUI::cb_dv_size(int xw)
{
    if( m_toplevel->resizable() )
    MxGUI::m_current->resize_canvas(xw, (2*xw)/3);
}

void MxGUI::cb_toggle(Fl_Menu_ *m, bool *flag)
{
    *flag = m->mvalue()->value()!=0;
    m_current->m_canvas->redraw();
}

void MxGUI::cb_save_view_to_file() { save_view_to_file(); }
void MxGUI::cb_load_view_from_file() { load_view_from_file(); }

bool MxGUI::save_view_to_file()
{
    fl_alert("This application has not defined a scheme for saving the view.");
    return false;
}

bool MxGUI::load_view_from_file()
{
    fl_alert("This application has not defined a scheme for loading the view.");
    return false;
}

int MxGUI::add_menu(const std::string& s, int key, Fl_Callback *cb, int flags)
{
    return m_menu_bar->add(s.c_str(), key, cb, this, flags);
}

int MxGUI::add_toggle_menu(const std::string& s, int key, bool& val, int flags)
{
    return m_menu_bar->add(s.c_str(), key, (Fl_Callback *)cb_toggle, &val,
             FL_MENU_TOGGLE|(val?FL_MENU_VALUE:0)|flags);
}

////////////////////////////////////////////////////////////////////////

#ifdef __CYGWIN__
extern "C"{
    extern int mainCRTStartup();
    int WinMainCRTStartup() { mainCRTStartup(); return 1; }
}
#endif

MxGUI::MxGUI()
{
    m_canvas = NULL;
    m_status_line = NULL;
    m_default_fps = 24.0f;
    m_target_fps = 0.0f;
    MxGUI::m_current = this;
}

Fl_Window *MxGUI::create_window(int xw, int yw, int pad)
{
    int yfill = 0;

    Fl_Window *w = new Fl_Window(xw+2*pad, 0);
    m_toplevel = w;
    w->box(FL_UP_BOX);

    m_menu_bar = new Fl_Menu_Bar(0, yfill, w->w(), 30);
    m_menu_bar->menu(m_menu_layout);
    yfill += m_menu_bar->h();

    add_upper_controls(yfill, pad);

    yfill += pad;
    m_canvas = new MxGLCanvas(pad, yfill, xw, yw);
    m_canvas->box(FL_DOWN_FRAME);
    m_canvas->attach_app(this);

    int glmode = 0;
    if(m_canvas->can_do(FL_RGB8))    glmode|=FL_RGB8;
    if(m_canvas->can_do(FL_DOUBLE))  glmode|=FL_DOUBLE;
    if(m_canvas->can_do(FL_DEPTH))   glmode|=FL_DEPTH;
    if(m_canvas->can_do(FL_ALPHA))   glmode|=FL_ALPHA;
    if(glmode)                       m_canvas->mode(glmode);

    yfill += m_canvas->h();

    add_lower_controls(yfill, pad);

    yfill += pad;
    m_status_line = new Fl_Output(pad, yfill, xw, 25);
    m_status_line->color(48);
    m_status_line->labeltype(FL_NO_LABEL);
    yfill += m_status_line->h();

    w->end();

    w->size(w->w(), yfill+pad);		// adjust window height
    w->resizable(*m_canvas);		// resize canvas with window

    // These are used by resize_canvas() to resize the window based on
    // the target size of the canvas.
    m_w_offset = w->w() - xw;
    m_h_offset = w->h() - yw;

    return w;
}

static
int arg_redirect(int argc, char **argv, int& index)
{
    MxGUI *app = MxGUI::m_current;
    return app?app->cmdline_option(argc, argv, index):0;
}

void MxGUI::initialize(int argc, char **argv, Fl_Menu_Item *m, int xw, int yw)
{
    Fl::visual(FL_RGB8);
    m_menu_layout = m?m:NULL;

    int index = 0;
    if( argv )
    Fl::args(argc, argv, index, arg_redirect);

    create_window(xw, yw);
    m_toplevel->label("Graphics Program");

    // Add dynamic entries
    typedef MxBinder<MxGUI> CB;
    std::string snap = "&File/Snapshot to/";
    std::string view = "&View/";
    std::string size = "&View/Display &size/";

    add_menu("&File/&New", FL_CTRL+'n', CB::to<&MxGUI::cb_new>);

#if defined(HAVE_LIBPNG)
    add_menu(snap+"&PNG",
        FL_CTRL+'p', CB::to_arg<&MxGUI::cb_snapshot, IMG_PNG>);
#endif
#if defined(HAVE_LIBTIFF)
    add_menu("&File/Snapshot to/&TIFF",
        FL_CTRL+'P', CB::to_arg<&MxGUI::cb_snapshot, IMG_TIFF>);
#endif
#if defined(HAVE_LIBJPEG)
    add_menu("&File/Snapshot to/&JPEG", 0, CB::to_arg<&MxGUI::cb_snapshot, IMG_JPEG>);
#endif
//    add_menu("&File/Snapshot to/PP&M", 0, CB::to_arg<&MxGUI::cb_snapshot, IMG_PNM>);

    add_menu("&File/E&xit", FL_CTRL+'q', CB::to<&MxGUI::cb_exit>);

    add_menu(view.append("Animation speed ..."), FL_CTRL+'r', CB::to<&MxGUI::cb_fps>);
    add_menu(view.append("&Animate"), FL_CTRL+'a', CB::to_menu<&MxGUI::cb_animate>, FL_MENU_TOGGLE);

    add_menu(view.append("Save view ..."), 0, CB::to<&MxGUI::cb_save_view_to_file>);
    add_menu(view.append("Load view ..."), 0, CB::to<&MxGUI::cb_load_view_from_file>);

    add_menu(size.append("320x240"),0, CB::to_arg<&MxGUI::cb_vga_size, 320>);
    add_menu(size.append("640x480"),0, CB::to_arg<&MxGUI::cb_vga_size, 640>);
    add_menu(size.append("800x600"),0, CB::to_arg<&MxGUI::cb_vga_size, 800>);
    add_menu(size.append("1024x768"),0,CB::to_arg<&MxGUI::cb_vga_size, 1024>,FL_MENU_DIVIDER);
    add_menu(size.append("720x480"),0, CB::to_arg<&MxGUI::cb_dv_size, 720>,FL_MENU_DIVIDER);

    add_menu(size.append("480x270"),0,   CB::to_arg<&MxGUI::cb_hdtv_size, 480>);
    add_menu(size.append("960x540"),0,   CB::to_arg<&MxGUI::cb_hdtv_size, 960>);
    add_menu(size.append("1920x1080"),0, CB::to_arg<&MxGUI::cb_hdtv_size, 1920>);

    if( argv )
    {
    if( index==argc )
        cmdline_file(NULL);
    else
        for(; index<argc; index++)
        cmdline_file(argv[index]);
    }
}

int MxGUI::run()
{
    m_toplevel->show();
    return Fl::run();
}

static
void cb_timeout(void *)
{
    MxGUI *app = MxGUI::m_current;

    if(!app || app->m_target_fps==0.0f)  return;

    app->update_animation();
    app->m_canvas->redraw();
    Fl::repeat_timeout(1/app->m_target_fps, cb_timeout);
}

void MxGUI::animate(bool will)
{
    if( will )
    {
    m_target_fps = m_default_fps;
    Fl::add_timeout(1/m_target_fps, cb_timeout);
    }
    else
    m_target_fps = 0;
}

int MxGUI::status(const char *fmt, ...)
{
    static char strbuf[1000];

    va_list args;
    va_start(args, fmt);
    int n = vsprintf(strbuf, fmt, args);

    m_status_line->value(strbuf);
    return n;
}

// bool MxGUI::snapshot_to_file(int format, const char *filenamep)
// {
//     canvas->make_current();
//     Fl::flush();

//     GLint vp[4]; glGetIntegerv(GL_VIEWPORT, vp);

//     glPushAttrib(GL_PIXEL_MODE_BIT);
//     glReadBuffer(GL_FRONT);
//     glPixelStorei(GL_PACK_ALIGNMENT, 1);

//     int nchan = 3;
//     if( format==IMG_PNG || format==IMG_TIFF )
// 	// Currently, only TIFF and PNG writers support the alpha channel
// 	nchan = 4;

//     ByteRaster img(vp[2]-vp[0], vp[3]-vp[1], nchan);
//     glReadPixels(vp[0],vp[1],vp[2],vp[3],
// 	    nchan==4 ? GL_RGBA : GL_RGB,
// 	    GL_UNSIGNED_BYTE,img.head());

//     glPopAttrib();
//     img.vflip();

//     if ( filenamep == NULL )
//     {
//         char msg[80], pat[8], name[16];
//         sprintf(msg, "Save %s Snapshot", image_type_name(format));
//         sprintf(pat, "*.%s", image_type_ext(format));
//         sprintf(name, "snap.%s", image_type_ext(format));

//         filenamep = fl_file_chooser(msg, pat, name);
//     }

//     if( filenamep )
// 	return write_image(filenamep, img, format);
//     else
// 	return false;
// }

void MxGUI::lock_size()
{
    m_toplevel->size_range(m_toplevel->w(), m_toplevel->h(),
            m_toplevel->w(), m_toplevel->h());
    m_toplevel->resizable(NULL);
}

void MxGUI::unlock_size()
{
    m_toplevel->resizable(*m_canvas);
    m_toplevel->size_range(100, 100, 0, 0);
}


void MxGUI::resize_canvas(int width, int height)
{
    m_toplevel->size(width+m_w_offset, height+m_h_offset);
    m_toplevel->redraw();
}

////////////////////////////////////////////////////////////////////////

void MxGUI::setup_for_drawing()
{
    cout << "MxGUI::setup_for_drawing()" << endl;
    glClearColor(0.65f, 0.65f, 0.65f, 0.0f);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1, 1, -1, 1, -1, 1);
}

void MxGUI::draw_contents()
{
    cout << "MxGUI::draw_contents()" << endl;
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glColor3f(0.0, 0.0, 0.0);
    glBegin(GL_LINES);
      glVertex2f(-1.0, 0.0);  glVertex2f(1.0, 0.0);
      glVertex2f(0.0, -1.0);  glVertex2f(0.0, 1.0);
    glEnd();
    
    glPopMatrix();
}

void MxGUI::update_animation()
{
}

bool MxGUI::mouse_down(int *where, int which)
{
    return false;
}

bool MxGUI::mouse_up(int *where, int which)
{
    return false;
}

bool MxGUI::mouse_drag(int *where, int *last, int which)
{
    return false;
}

bool MxGUI::key_press(int key)
{
    return false;
}

int MxGUI::cmdline_option(int argc, char **argv, int& index)
{
    return 0;
}

void MxGUI::cmdline_file(const char *file)
{
}

#endif /* HAVE_FLTK */
