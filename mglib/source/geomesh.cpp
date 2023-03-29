#include "auxiliar.hpp"
#include "geomesh.hpp"
#include "general_math.hpp"
#include "path_tracer.hpp"
#include "geom3d.hpp"

#include <float.h>

bool DEBUG = false;

GeoTriMesh::GeoTriMesh()
: TriMesh(), ERROR_TOLERANCE(0.001)
{
    initialize();
}

GeoTriMesh::GeoTriMesh(TriMesh *m)
: TriMesh(*m), ERROR_TOLERANCE(0.001)
{
    initialize();
}

GeoTriMesh::~GeoTriMesh()
{
    reset_distance();

    if (graph) delete graph;
}
 
void GeoTriMesh::initialize()
{
    graph = new ManifoldGraphT(m_vertex.size(), m_face);
    prop_step_size = 1;
    num_prop_edge = 0;
    is_first_propagation = true;
    stop_prop_number = 283;

    geo_distance.resize(m_vertex.size(), DBL_MAX);
}

void GeoTriMesh::reset_distance()
{
    is_first_propagation = true;
    num_prop_edge = 0;

    clear_edge_map();
    i_queue.reset();
    geo_distance.clear();
    geo_distance.resize(m_vertex.size(), DBL_MAX);
    last_step.reset();
}

void GeoTriMesh::clear_edge_map()
{
    EdgeStructMap::iterator esit = edge_map.begin();
    for ( ; esit != edge_map.end(); esit++) {
        EdgeStruct& es = (*esit).second;
        EdgeStruct::iterator eit = es.begin();
        for ( ; eit != es.end(); eit++) {
            delete (*eit); (*eit) = NULL;
        }
        es.clear();
    }
    edge_map.clear();
}

void GeoTriMesh::compute_geodesic(int selected_v, ScalarVector& distance, KnotVectorVector& path)
{
    reset_distance();

    if (is_first_propagation)
        init_propagation(selected_v);

    // propagate intervals and compute distance field
    while (propagate_once());
    // trace the geodesic paths
    backtracing();

    // output the results
    distance.resize(m_vertex.size());
    copy(geo_distance.begin(), geo_distance.end(), distance.begin());

    path.resize(m_vertex.size());
    copy(geo_path.begin(), geo_path.end(), path.begin());
}

bool GeoTriMesh::init_propagation(int selected_v)
{
    Vec3& source_v = m_vertex[selected_v];

    // initialize the first elements in the priority queue
    HandleSet f_edges, a_edges;
    graph->collect_vertex_front_edge(selected_v, f_edges);
    graph->collect_vertex_adj_edge_with_sym(selected_v, a_edges);

    HandleSet::iterator hit;
    
    // initialize intervals for the adjacent edges and insert them into edge map
    for (hit = a_edges.begin(); hit != a_edges.end(); hit++)
    {
        Interval* iv = create_init_interval(*hit, source_v);
        edge_map[iv->handle()].push_back(iv);
    }

    // initialize interval for the front edges, 
    // inser them into edge map, and insert them into priority queue
    for (hit = f_edges.begin(); hit != f_edges.end(); hit++)
    {
        Interval* iv = create_init_interval(*hit, source_v);
        edge_map[iv->handle()].push_back(iv);
        i_queue.insert(iv);
    }

    geo_distance[m_selected_vertex] = 0.f;
    for (EdgeStructMap::iterator eit = edge_map.begin(); eit != edge_map.end(); eit++)
        update_sfield((*eit).first, (*eit).second);

    is_first_propagation = false;
    return true;
}

// propagate one interval
bool GeoTriMesh::propagate_once()
{
    if (DEBUG)
        DEBUG = DEBUG;

    apply_last_step();

    if (i_queue.empty()) return false;

    Interval* from_iv = (Interval*) i_queue.extract();

    // propagate interval
    EdgeStructMap new_ivs;
    propagate_interval(from_iv, new_ivs);

    // merge interval
    for (EdgeStructMap::iterator nit = new_ivs.begin(); nit != new_ivs.end(); nit++)
    {
        Handle to_e = (*nit).first;
        EdgeStruct& n_ivs = (*nit).second;
        EdgeStruct& o_ivs = edge_map[to_e];
        EdgeStruct copy_o_ivs;
        copy_vector(o_ivs, copy_o_ivs);

        if (o_ivs.empty())
        {
            copy_o_ivs.insert_intervals(n_ivs);

            // store the modified intervals
            last_step.modified_ivs.push_back(make_pair(from_iv, Propagation(copy_o_ivs)));
        }
        else 
        {
            for (EdgeStruct::iterator nnit = n_ivs.begin(); nnit != n_ivs.end(); nnit++)
            {
                Propagation propag;

                // merge intervals
                merge_interval(*nnit, copy_o_ivs, propag);

                // if (gui.will_approximate_interval)
                //     approximate_interval(propag);

                copy_o_ivs.replace_interval_list(propag.altered);
                // store the modified intervals
                last_step.modified_ivs.push_back(make_pair(from_iv, propag));
            }
        }
    }
    from_iv->set_finished();
    return true;
}

void GeoTriMesh::apply_last_step()
{
    if (last_step.cur_iv == NULL) return;

    vector< pair<Interval*, Propagation> >::iterator pit = last_step.modified_ivs.begin();
    for ( ; pit != last_step.modified_ivs.end(); pit++) 
    {
        Interval* from_iv = (*pit).first;
        Propagation& propag = (*pit).second;
        Handle to_e = propag.altered.front()->handle();

        EdgeStruct& o_ivs = edge_map[to_e];
        o_ivs.replace_interval_list(propag.altered);

        update_heap(from_iv, propag);
        update_sfield(to_e, propag.altered);
    }
    last_step.reset();
}

void GeoTriMesh::backtracing()
{
    PathTracer tracer(*this, geo_path);
    tracer.trace_path();
}

Vec2 GeoTriMesh::project_p(Handle e, int to_face, const Vec3& p) const 
{
    Vec2 proj_p;
    Vec3 x_axis = edge_vector(e);
    Vec3 p_vec = p - m_vertex[e->Org()];
    Vec3 y_axis = m_fnormal[to_face] ^ x_axis;
    unitize(y_axis);

    proj_p[0] = p_vec*x_axis;
    proj_p[1] = p_vec*y_axis;

    return proj_p;
}

Vec3 GeoTriMesh::unproject_p(Handle e, int to_face, const Vec2& s_xy) const
{
    Vec3 s_xyz;
    Vec3 org = m_vertex[e->Org()];

    Vec3 x_axis = edge_vector(e);
    Vec3 y_axis = m_fnormal[to_face] ^ x_axis;
    unitize(y_axis);
    
    s_xyz = s_xy[0]*x_axis + s_xy[1]*y_axis;
    return org + s_xyz;
}

bool GeoTriMesh::is_iv_b0_end(bool tau, Interval* iv, double e_length) 
{
    return (!tau) ? FEQ(iv->get_b0(), 0) : FEQ(iv->get_b1(), e_length);
}

bool GeoTriMesh::is_iv_b1_end(bool tau, Interval* iv, double e_length) 
{
    return (!tau) ? FEQ(iv->get_b1(), e_length) : FEQ(iv->get_b0(), 0);
}

void GeoTriMesh::update_sfield(Handle to_e, EdgeStruct& ivs)
{
    if (to_e != NULL && !ivs.empty())
    {
        int v_endp1 = to_e->Org();
        int v_endp2 = to_e->Dest();
        double elen = edge_length(to_e);

        // for DEBUG
        double temp_smallest_b0 = ivs.get_smallest_b0(), temp_largest_b1 = ivs.get_largest_b1();
        if (FEQ(0.f, ivs.get_smallest_b0()))
        {
            double old_org = geo_distance[v_endp1];
            double new_org = ivs.get_distance_end0();
            if (FLT(new_org, old_org)) geo_distance[v_endp1] = new_org;
        }
        if (FEQ(elen, ivs.get_largest_b1()))
        {
            double old_dest = geo_distance[v_endp2];
            double new_dest = ivs.get_distance_end1();
            if (FLT(new_dest, old_dest)) geo_distance[v_endp2] = new_dest;
        }
    }
}

Interval* GeoTriMesh::create_init_interval(Handle e, Vec3& source)
{
    Handle eh;
    double b0 = 0.0;
    double b1 = edge_length(e);
    bool tau;

    if (e->Org() > e->Dest() && e->Sym() != NULL)
    { eh = e->Sym(); tau = true; }
    else 
    { eh = e; tau = false; }

    double d0 = norm(m_vertex[eh->Org()]-source);
    double d1 = norm(m_vertex[eh->Dest()]-source);
    Interval* iv = new Interval(eh, b0, b1, d0, d1, 0.f, tau);

    return iv;
}

Interval* GeoTriMesh::create_interval(Handle e, double sigma, Vec3& source, double e_l)
{
    double e_length = (e_l<0) ? edge_length(e) : e_l;
    double d0 = norm(source-m_vertex[e->Org()]);
    double d1 = norm(source-m_vertex[e->Dest()]);

    return create_interval(e, 0, e_length, d0, d1, sigma);
}

Interval* GeoTriMesh::create_interval(Handle e, double sigma, double b0, double b1, Vec3& source, bool inverted, double e_l)
{
    double d0 = norm(source-edge_point(e, b0));
    double d1 = norm(source-edge_point(e, b1));

    return create_interval(e, b0, b1, d0, d1, sigma, inverted);
}

Interval* GeoTriMesh::create_interval(Handle e, double b0, double b1, double d0, double d1, double sigma, bool inverted)
{
    Handle eh = e; 
    Interval* iv;
    
    if (e->Org() <= e->Dest()) {
        // not inverted
        iv = new Interval(eh, b0, b1, d0, d1, sigma, inverted);
    } else {
        // invert the interval attributes
        eh = e->Sym();
        double e_length = edge_length(eh);
        double bb0 = e_length - b1;
        double bb1 = e_length - b0;
        iv = new Interval(eh, bb0, bb1, d1, d0, sigma, true);
    }
    return iv;
}

Vec2 GeoTriMesh::get_s_coord(const Interval* iv) const
{
    Vec2 s;
    double elength = edge_length(iv->handle());
    double b0, b1, r0, r1;
    if (!iv->get_tau())
    {
        b0 = iv->get_b0(); b1 = iv->get_b1();
        r0 = iv->get_d0(); r1 = iv->get_d1();
    } else {
        b0 = elength - iv->get_b1(); b1 = elength - iv->get_b0();
        r0 = iv->get_d1(); r1 = iv->get_d0();	
    }

    double d = b1-b0;
    double temp = d*d-r1*r1+r0*r0;
    double y_square = (4.0*d*d*r0*r0 - temp*temp)/(4.0*d*d);

    while (FLT(y_square, 0.0)) {
        if (r0 < r1)
            r0 += (r0*1e-6); 
        else 
            r1 += (r1*1e-6);
        temp = d*d-r1*r1+r0*r0;
        y_square = (4.0*d*d*r0*r0 - temp*temp)/(4.0*d*d);
    }

    s[0] = b0 + temp/(2.0*d);
    s[1] = FEQ(y_square, 0.0) ? 0.0 : sqrt(y_square);

    return s;
}

Vec2 GeoTriMesh::get_s_coord_straight(const Interval* iv) const
{
    Vec2 s;
    double elength = edge_length(iv->handle());
    double b0, b1, r0, r1;
    b0 = iv->get_b0(); b1 = iv->get_b1();
    r0 = iv->get_d0(); r1 = iv->get_d1();

    double d = b1-b0;
    double temp = d*d-r1*r1+r0*r0;
    double y_square = (4.0*d*d*r0*r0 - temp*temp)/(4.0*d*d);

    while (FLT(y_square, 0.0)) {
        if (r0 < r1)
            r0 += (r0*1e-6); 
        else 
            r1 += (r1*1e-6);
        temp = d*d-r1*r1+r0*r0;
        y_square = (4.0*d*d*r0*r0 - temp*temp)/(4.0*d*d);
    }

    s[0] = b0 + temp/(2.0*d);
    s[1] = FEQ(y_square, 0.0) ? 0.0 : sqrt(y_square);

    return s;
}

double GeoTriMesh::edge_length(const Handle& e) const
{
    return norm(m_vertex[e->Org()] - m_vertex[e->Dest()]); 
}

Vec3 GeoTriMesh::edge_vector(const Handle& e) const
{
    Vec3 v = m_vertex[e->Dest()] - m_vertex[e->Org()];
    unitize(v);
    return v;
}

Vec3 GeoTriMesh::edge_point(const Handle& handle, double interv) const 
{
    Vec3 p_org = m_vertex[handle->Org()];
    Vec3 p_dest = m_vertex[handle->Dest()];
    Vec3 p_vec = p_dest - p_org;
    double ratio = 0, e_length = norm(p_vec);
    if (e_length > 0)
        ratio = interv / e_length;

    return p_org + ratio * p_vec;
}

Vec3 GeoTriMesh::edge_point_unit(const Handle& handle, double ratio) const 
{
    Vec3 p_org = m_vertex[handle->Org()];
    Vec3 p_dest = m_vertex[handle->Dest()];
    Vec3 p_vec = p_dest - p_org;

    return p_org + ratio * p_vec;
}

double GeoTriMesh::edge_interval(const Vec2& p0, const Vec2& p1, const Vec2& point) const
{
    if (GMath::are_ordered_along_line(point, p0, p1))
        return -norm(p0-point);
    else 
        return norm(p0-point);
}

double GeoTriMesh::edge_interval_unit(const Vec2& p0, const Vec2& p1, const Vec2& point) const
{
    double e_length = norm(p1-p0);
    if (GMath::are_ordered_along_line(point, p0, p1))
        return -norm(p0-point)/e_length;
    else 
        return norm(p0-point)/e_length;
}


double GeoTriMesh::distance_at(const Interval* iv, const double inter) const
{
    Handle e = iv->handle();
    double dist = dvalue_at(iv, inter);

    return iv->sigma + dist;
}

double GeoTriMesh::dvalue_at(const Interval* iv, const double inter) const
{
    Handle e = iv->handle();
    double dist, e_length = edge_length(e);
    Vec2 s = get_s_coord_straight(iv);
    dist = norm(s-Vec2(inter, 0));

    return dist;
}

void GeoTriMesh::update_heap(Interval* from_iv, Propagation& propag)
{
    // delete the discarted intervals from the heap
    IntervalList::iterator lit;
    for (lit = propag.deleted.begin(); lit != propag.deleted.end(); lit++)
    {
        i_queue.remove(*lit); delete (*lit); (*lit) = NULL;
    }

    // update the altered intervals in the heap
    for (lit = propag.altered.begin(); lit != propag.altered.end(); lit++)
    {
        if ((*lit) != from_iv)
            if ((*lit)->is_in_heap()) 
                i_queue.update(*lit);
            else if (!(*lit)->is_finished()) 
                i_queue.insert(*lit);
    }
}

void GeoTriMesh::update_heap(Interval* from_iv, IntervalList& altered)
{
    IntervalList::iterator lit;
    // update the altered intervals in the heap
    for (lit = altered.begin(); lit != altered.end(); lit++)
        if ((*lit) != from_iv)
            if ((*lit)->is_in_heap()) 
                i_queue.update(*lit);
            else if (!(*lit)->is_finished()) 
                i_queue.insert(*lit);
}

void GeoTriMesh::set_error_threshold(double thres)
{
    ERROR_TOLERANCE = thres;
}

Vec3 GeoTriMesh::vertex_normal(int vid)
{
    FaceSet fring;
    graph->collect_face_ring(vid, fring);
    Vec3 vn;
    FaceSet::iterator fit = fring.begin();
    for ( ; fit != fring.end(); fit++) {
        cout << "fid: " << *fit << endl;
        vn += m_fnormal[*fit];
    }

    unitize(vn);
    return vn;
}