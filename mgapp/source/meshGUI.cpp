
/************************************************************************
 * File description: Main file for the mesh geodesy GUI related functions
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/


#include "mxGUI.hpp"
#include "meshGUI.hpp"
#include "vec4.hpp"
#include "gltools.hpp"

#include <FL/Fl_File_Chooser.H>

MeshGUI::MeshGUI() : MxGUI(), 
        m_draw_mode(Draw_mode_solid), 
        m_selection_mode(Noselect),
        m_selected_vertex(-1),
        m_will_draw_vertices(false),
        m_will_draw_bbox(false), 
        m_will_draw_surface_fnormal(true), 
        m_will_draw_mesh(false),
        m_will_draw_geodesic_distance(false),
        m_will_draw_geodesic_path(false),
        m_grid_period(0x10),
        m_mesh(nullptr),
        m_mls(nullptr),
        m_displayed_face_id(0),
        m_display_face_increment(400)
{
    m_mat = new Material();
    m_obj = gluNewQuadric();

    default_texture();
}

MeshGUI::~MeshGUI() 
{
    if (m_mesh)     delete m_mesh;
    if (m_mls)      delete m_mls;

    if (m_texture)  delete m_texture;
    if (m_mat)      delete m_mat;
    if (m_obj)      gluDeleteQuadric(m_obj);
}

void MeshGUI::initialize(int argc, char* argv[])
{
    LOG(INFO) << "MeshGUI::initialize";
    MxGUI::initialize(argc, argv);

    // Set up the menu bar
    add_menu_item("&File/Open", FL_CTRL+'o', (Fl_Callback *)cb_open_file, NULL);
    add_menu_item("&File/Open", FL_CTRL+'o', (Fl_Callback *)cb_open_file, NULL);
    add_menu_item("&File/Save", FL_CTRL+'s', (Fl_Callback *)cb_save_file, NULL);

    add_toggle_menu("&Draw/Vertices",              0, m_will_draw_vertices);
    add_toggle_menu("&Draw/Surface + face normal", 0, m_will_draw_surface_fnormal);
    add_toggle_menu("&Draw/Mesh",                  0, m_will_draw_mesh);
    add_toggle_menu("&Draw/Bounding box", FL_CTRL+'b', m_will_draw_bbox);
    add_toggle_menu("&Draw/Geodesic distance",     0, m_will_draw_geodesic_distance);
    add_toggle_menu("&Draw/Geodesic path",         0, m_will_draw_geodesic_path);

    add_menu_item("&Tools/Default texture",        0, (Fl_Callback *)cb_load_default_texture, NULL);
    add_menu_item("&Tools/Load texture",           0, (Fl_Callback *)cb_load_texture, NULL);
    add_menu_item("&Tools/Save geodesic distance", 0, (Fl_Callback *)cb_save_distance, NULL);
    add_menu_item("&Tools/Load geodesic distance", 0, (Fl_Callback *)cb_load_distance, NULL);
}

int MeshGUI::add_menu_item(const char* name, int key, Fl_Callback *f, void* val, int flags) 
{
    return MxGUI::m_menu_bar->add(name, key, f, val);
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
        if (gui.m_mesh)
            gui.m_mesh->write_to_file(output_filename);
    }
}

void MeshGUI::cb_load_default_texture()
{
    gui.default_texture();
    gui.setup_texture();
    gui.m_canvas->redraw();
}

void MeshGUI::cb_load_texture()
{
    const char *filename = fl_file_chooser("Load texture:", "*.png", "");
    if (!filename) return;

    ByteRaster *img = read_image(filename);
    if( img )
    {
        if (gui.m_texture) delete gui.m_texture;
        gui.m_texture = img;
        gui.m_canvas->make_current();
        gui.setup_texture();
    }
}

void MeshGUI::cb_save_distance() 
{
    string input_filename = fl_file_chooser("Select input file", "*{.dist}", NULL);

    if (!input_filename.empty())
        gui.m_mls->save_distances(input_filename);
}

void MeshGUI::cb_load_distance() 
{
    string output_filename = fl_file_chooser("Select output file", "*{.dist}", NULL);

    if (!output_filename.empty())
        gui.m_mls->load_distances(output_filename);
}

void MeshGUI::load_mesh(const string& filename)
{
    LOG(INFO) << "MeshGUI::load_mesh " << filename;
    if (m_mls) delete m_mls;
    TriMesh *tri = new TriMesh();

    tri->read_from_file(filename);
    tri->initialize();
    tri->compute_bbox(m_bb_min, m_bb_max);
    m_mesh = new GeoTriMesh(tri);

    m_mls = new MLS(m_mesh);

    reset_camera();
    m_canvas->redraw();
}

void MeshGUI::setup_for_drawing() 
{
    LOG(INFO) << "MeshGUI::setup_for_drawing";
    glClearColor(1.f, 1.f, 1.f, 0.0f); // white
    // glClearColor(26.f/255.f, 34.f/255.f, 40.f/255.f, 0.0f);  // dark blue
    
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

    m_mat->load_standard();

    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    setup_texture();
}

void MeshGUI::setup_texture()
{
    LOG(INFO) << "MeshGUI::setup_texture";
    GLenum fmt;
    switch( m_texture->channels() )
    {
    case 1:  fmt=GL_LUMINANCE; break;
    case 3:  fmt=GL_RGB; break;
    case 4:  fmt=GL_RGBA; break;
    default:
         cerr << "Sorry, but I need a valid texture!" << endl;
         exit(1);
    }

    // setting up the texture parameters
    const GLenum TEX = GL_TEXTURE_2D;

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

#if 1
    glTexImage2D(GL_TEXTURE_2D, 0, m_texture->channels(), m_texture->width(), m_texture->height(), 0,
        fmt, GL_UNSIGNED_BYTE, m_texture->head());
#else
    gluBuild2DMipmaps(TEX, m_texture->channels(), m_texture->width(), m_texture->height(), 
        fmt, GL_UNSIGNED_BYTE, m_texture->head());
#endif

    glTexParameterf(TEX, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(TEX, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(TEX, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(TEX, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

void MeshGUI::default_texture()
{
    LOG(INFO) << "MeshGUI::default_texture";
    if( m_texture )  delete m_texture;

    // Set up an initial bogus texture
    m_texture = new ByteRaster(512, 512, 1);
    const int block = m_grid_period;
    for(int i=0; i<m_texture->width(); ++i)
        for(int j=0; j<m_texture->height(); ++j)
        {
            int g = ( ((i&block)==0) ^ ((j&block)==0) )*255;
            m_texture->pixel(i,j)[0] = g;
        }
    }

void MeshGUI::draw_contents()
{
    LOG(INFO) << "MeshGUI::draw_contents";
    begin_redraw();

    glPushMatrix();
    default_redraw();
    glPopMatrix();
    
    end_redraw();
}

void MeshGUI::begin_redraw() 
{
    LOG(INFO) << "MeshGUI::begin_redraw";
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    
    if( m_selection_mode == Noselect ) {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
    }
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    apply_camera();
    
    m_ball.apply_transform();
}

void MeshGUI::default_redraw() 
{
    LOG(INFO) << "MeshGUI::default_redraw";
    if (m_selected_vertex != -1)     draw_selection();
    if (m_will_draw_bbox)            draw_bbox();
    if (m_will_draw_surface_fnormal) draw_surface_fnormal();
    if (m_will_draw_mesh)            draw_mesh();
    if (m_will_draw_vertices)        draw_vertices();

    if (m_will_draw_geodesic_distance)	draw_geodesic_distance();
    if (m_will_draw_geodesic_path)		draw_geodesic_path();
}

void MeshGUI::end_redraw() 
{
    LOG(INFO) << "MeshGUI::end_redraw";
    m_ball.unapply_transform();
}

void MeshGUI::apply_camera() 
{
    float aspect = (float)m_canvas->w() / (float)m_canvas->h();
    camera_lookat(m_bb_min, m_bb_max, aspect);
}

void MeshGUI::reset_camera() 
{
    if (m_mesh) m_mesh->compute_bbox(m_bb_min, m_bb_max);
    Vec3 ctr = (m_bb_max + m_bb_min)/2.0;
    double radius = norm(m_bb_max - ctr);
    m_ball.bounding_sphere(ctr, radius);
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

void MeshGUI::draw_bbox() 
{
    LOG(INFO) << "MeshGUI::draw_bbox";
    glPushAttrib(GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(0.f, 0.f, 0.f);
    draw_box(m_bb_min, m_bb_max);
    glPopAttrib();
}

void MeshGUI::draw_box(const Vec3f& min, const Vec3f& max) 
{
    LOG(INFO) << "MeshGUI::draw_box";
    glBegin(GL_LINE_LOOP);
    glVertex3f(min[0],min[1],min[2]); glVertex3f(min[0],max[1],min[2]);
    glVertex3f(max[0],max[1],min[2]); glVertex3f(max[0],min[1],min[2]);
    glEnd();
    
    glBegin(GL_LINE_LOOP);
    glVertex3f(min[0],min[1],max[2]); glVertex3f(min[0],max[1],max[2]);
    glVertex3f(max[0],max[1],max[2]); glVertex3f(max[0],min[1],max[2]);
    glEnd();
    
    glBegin(GL_LINES);
    glVertex3f(min[0],min[1],min[2]); glVertex3f(min[0],min[1],max[2]);
    glVertex3f(min[0],max[1],min[2]); glVertex3f(min[0],max[1],max[2]);
    glVertex3f(max[0],max[1],min[2]); glVertex3f(max[0],max[1],max[2]);
    glVertex3f(max[0],min[1],min[2]); glVertex3f(max[0],min[1],max[2]);
    glEnd();
}

void MeshGUI::draw_vertices() 
{
    LOG(INFO) << "MeshGUI::draw_vertices";
    glVertexPointer(3, GL_DOUBLE, 0, &m_mesh->m_vertex[0]);

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(0.f, 0.f, 0.f);

    int psize = m_mesh->m_vertex.size();
    glBegin(GL_POINTS);
    for (int i=0; i<psize; i++) {
        glArrayElement(i);
    }
    glEnd();
    
    glPopAttrib();    
}

void MeshGUI::draw_mesh() 
{
    LOG(INFO) << "MeshGUI::draw_mesh";
    glDisable(GL_POLYGON_OFFSET_FILL);
    
    glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
    glVertexPointer(3, GL_DOUBLE, 0, &m_mesh->m_vertex[0]);
    
    glDisable(GL_LIGHTING);
    glColor3f(0.f, 0.f, 0.f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    TriMesh::FaceList& face = m_mesh->m_face;
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
    
     if (!m_will_draw_surface_fnormal) 
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
    LOG(INFO) << "MeshGUI::draw_surface_fnormal";
    glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
    glVertexPointer(3, GL_DOUBLE, 0, &m_mesh->m_vertex[0]);
    
    if (m_will_draw_mesh) 
    {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
    }
    
    if (m_mesh->m_fnormal.empty())
        m_mesh->compute_fnormal();

    TriMesh::FaceList& face = m_mesh->m_face;
    int fsize = (m_target_fps == 0) ? face.size() : m_displayed_face_id;
    if (m_draw_mode == Draw_mode_solid) 
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
    else if (m_draw_mode == Draw_mode_wireframe) 
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

void MeshGUI::draw_selection()
{
    LOG(INFO) << "MeshGUI::draw_selection";
    const float sball_radius = 0.01f;

    glPushAttrib(GL_ENABLE_BIT);

    m_mat->load_red();
    if (m_selected_vertex != -1)
    {
        Vec3& v = m_mesh->m_vertex[m_selected_vertex];
        // draw a ball in the vertex position
        glPushMatrix();
        glTranslated(v[0], v[1], v[2]);
        gluSphere(m_obj, sball_radius, 10, 10);
        glPopMatrix();
    }
    m_mat->load_standard();

    glPopAttrib();
}

void MeshGUI::draw_geodesic_distance()
{
    LOG(INFO) << "MeshGUI::draw_geodesic_distance";
    if (m_mls->distances.empty()) return;

    glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
    glVertexPointer(3, GL_DOUBLE, 0, &m_mesh->m_vertex[0]);

    const TriMesh::FaceList& face = m_mesh->m_face;

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glEnable(GL_DEPTH_TEST);	// Enable the Z-Buffer
    glEnable(GL_TEXTURE_2D);

    int size = face.size();
    GeoTriMesh::ScalarVector& dist = m_mls->distances;
    for(int i=0; i<size; i++)
    {
        const Face& f = face[i];
        glBegin(GL_TRIANGLES);
        setup_face_state(i);
        glTexCoord2f(dist[f[0]], 1.f);
        glArrayElement(f[0]);
        glTexCoord2f(dist[f[1]], 1.f);
        glArrayElement(f[1]);
        glTexCoord2f(dist[f[2]], 1.f);
        glArrayElement(f[2]);
        glEnd();
    }

    glPopAttrib();
}

void MeshGUI::draw_geodesic_path()
{
    LOG(INFO) << "MeshGUI::draw_geodesic_path";
    if (m_mls->paths.empty()) return;

    glPushAttrib(GL_ENABLE_BIT|GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(0.3f,0.4f,0.5f);
    glLineWidth(0.5f);

    for (auto& path : m_mls->paths ) {
        MLS::KnotVector::iterator vit = path.begin();
        glBegin(GL_LINE_STRIP);
        for ( ; vit != path.end(); vit++)
            glVertex3dv(m_mesh->edge_point((*vit).first, (*vit).second));
        glEnd();
    }
    glPopAttrib();
}

void MeshGUI::setup_face_state(int fid) 
{
    if (fid < m_mesh->m_fnormal.size())
        glNormal3dv(m_mesh->m_fnormal[fid]);
}

bool MeshGUI::mouse_down(int *where, int which) 
{
    LOG(INFO) << "MeshGUI::mouse_down";
    if (which == 1 && Fl::event_state(FL_SHIFT)) 
    {
        int old = m_selected_vertex;
        m_selected_vertex = pick_vertex(where);
        cout << "selected vertex: " << m_selected_vertex << endl;
        m_canvas->redraw();
        return old != m_selected_vertex;
    } 
    return m_ball.mouse_down(where, which);
}

bool MeshGUI::mouse_drag(int *where, int *last, int which) 
{
    if (which == 1 && Fl::event_state(FL_ALT))
        return m_ball.mouse_drag(where, last, 3);
    else if (which == 1 && Fl::event_state(FL_SHIFT))
        return m_ball.mouse_drag(where, last, 2);
    return m_ball.mouse_drag(where, last, which);
}

bool MeshGUI::mouse_up(int *where, int which) 
{
    return m_ball.mouse_up(where, which);
}

void MeshGUI::animate(bool will)
{
    MxGUI::animate(will);
    m_target_fps = 1;
}

void MeshGUI::update_animation()
{
    m_displayed_face_id += m_display_face_increment;
    m_displayed_face_id %= m_mesh->m_face.size();

    int where[2] = {350, 282};
    int last[2] = {348,282};
    m_ball.mouse_drag(where, last, 1);
    apply_camera();
    m_canvas->redraw();
}

bool MeshGUI::key_press(int key) 
{    
    LOG(INFO) << "MeshGUI::key_press";
    switch (key) {
    case FL_Up:
        up_frequency(); 
        m_canvas->redraw(); 
        break;
    case FL_Down:
        down_frequency(); 
        m_canvas->redraw(); 
        break;
    case 'w':
        m_draw_mode = Draw_mode_wireframe;
        m_canvas->redraw();
        break;
    case 's':
        m_draw_mode = Draw_mode_solid;
        m_canvas->redraw();
        break;
    case 'n':
        m_mesh->normalize();
        reset_camera();
        m_canvas->redraw();
        break;
    case 'd':
        m_selected_vertex = (m_selected_vertex == -1) ? 0 : m_selected_vertex;
        m_mls->compute_distances(m_selected_vertex);
        cout << "Calculating geodesic distances" << endl;
        m_canvas->redraw();
        break;
    case 'p':
        m_mls->sort_faces_by_distance();
        m_canvas->redraw();
        break;
    }
    return true;
}

int MeshGUI::pick_vertex(int where[2]) 
{
    LOG(INFO) << "MeshGUI::pick_vertex";
    GLuint buffer[128];
    double radius = 16.0;
    
    m_selection_mode = Vselect;
    m_canvas->make_current();
    
    begin_opengl_pick(where, radius, buffer, 128);
    begin_redraw();
    draw_for_selection();
    end_redraw();
    m_selection_mode = Noselect;
    
    return complete_opengl_pick(buffer);    
}

void MeshGUI::draw_for_selection() 
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_DOUBLE, 0, &m_mesh->m_vertex[0]);

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
    glDisable(GL_LIGHTING);
    glColor3f(0.f, 0.f, 0.f);
    glPointSize(5.f);
    
    if (m_selection_mode == Vselect) {
        int psize = m_mesh->m_vertex.size();
        for (int i=0; i<psize; i++) {
            glLoadName(i);
            glBegin(GL_POINTS);
            glArrayElement(i);
            glEnd();
        }
    }
    glEnable(GL_LIGHTING);
    glPopAttrib();    
}

void MeshGUI::up_frequency()
{
    m_grid_period /= 2;
    default_texture();
    m_canvas->make_current();
    setup_texture();
}

void MeshGUI::down_frequency()
{
    m_grid_period *= 2;
    default_texture();  
    m_canvas->make_current();
    setup_texture();
}