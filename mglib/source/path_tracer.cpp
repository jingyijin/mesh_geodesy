/************************************************************************
 * File description: Path tracer class for the path tracing algorithm
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/


#include "path_tracer.hpp"

PathTracer::PathTracer(GeoTriMesh& m, PathKnotVectorVector& k_path) 
: m_mesh(m), m_path_knot(k_path), m_vindex(0), m_state(Init_state)
{}

void PathTracer::clear_path()
{
    vector<PathKnotVector>::iterator kit = m_path_knot.begin();
    for ( ; kit != m_path_knot.end(); kit++)
        (*kit).clear();
    m_path_knot.clear();
}

void PathTracer::trace_path()
{
    clear_path();

    int vsize = m_mesh.m_vertex.size();
    m_path_knot.resize(vsize);
    m_vindex = 0;

    while (m_vindex < vsize)
        trace_once();
}

void PathTracer::trace_once()
{
    switch (m_state) {
    case Init_state:
    {
        m_target_vertex = m_vindex;
        // store the path knot and coordinates
        store_path_knot(m_target_vertex, m_target_vertex);
        m_next_vertex = m_target_vertex;
        m_prev_e = NULL;
        if (m_target_vertex == m_mesh.m_selected_vertex) 
        {
            m_state = Complete_state;
            break;
        }
        else 
            m_state = Tracing_ring_build_state;
    }
    case Tracing_ring_build_state:
    {
        build_current_ring(m_next_vertex, m_prev_e, m_current_ring);
        m_state = Tracing_ring_trace_state;
    }
    case Tracing_ring_trace_state:
    {
        m_min_iv = get_one_interval(m_current_ring, m_start_coord);
        if (m_min_iv == NULL) 
            m_state = Complete_state;
        else {
            m_prev_e = m_min_e;
            m_min_e = m_min_iv->handle();

            int new_vertex = get_next_vertex(m_min_iv);
            if (new_vertex >= 0 && new_vertex != m_next_vertex) 
            {
                // store the path knot and coordinates
                store_path_knot(m_target_vertex, new_vertex, m_min_e);

                m_next_vertex = new_vertex;
                m_current_ring.clear();
                if (FEQ(m_min_iv->get_sigma(), 0.0))
                    m_state = Complete_state;
                else 
                    m_state = Tracing_ring_build_state;
            } 
            else if (new_vertex == m_next_vertex)
            {
                flag_current_ring(m_current_ring, m_min_e);
                m_state = Tracing_ring_trace_state;
            }
            else 
            {
                flag_current_ring(m_current_ring, m_min_e);
                if (!m_min_iv->get_tau()) 
                    m_cur_e = m_min_e;
                else
                    m_cur_e = m_min_e->Sym();
                m_to_face = m_cur_e->Lface();
                m_s = m_mesh.get_s_coord(m_min_iv);
                m_ray = Ray<Vec2>(m_start_coord, m_s);

                Interval* next_iv = trace_ray_once(m_to_face, m_cur_e, m_ray, m_inter);
                if (next_iv != NULL) 
                {
                    Handle next_e = next_iv->handle();

                    if ((m_current_ring.find(next_e) != m_current_ring.end() || 
                         m_current_ring.find(next_e->Sym()) != m_current_ring.end()))
                    {
                        m_state = Tracing_ring_trace_state;
                    }
                    else 
                    {
                        // exiting the tracing in a ring 
                        m_current_ring.clear();

                        // store the path knot and coordinate
                        store_path_knot(m_target_vertex, m_inter.e, m_inter.inter);

                        // prepare for the next step tracing
                        m_min_iv = next_iv;
                        m_prev_e = m_min_e;
                        m_min_e = m_min_iv->handle();
                        if (m_min_iv->handle()->Org() == m_inter.e->Org())
                        {
                            if (!m_min_iv->get_tau())
                                m_start_coord = Vec2(m_inter.inter, 0);
                            else 
                                m_start_coord = Vec2(m_mesh.edge_length(m_inter.e)-m_inter.inter, 0);
                        } 
                        else 
                        {
                            if (!m_min_iv->get_tau())
                                m_start_coord = Vec2(m_mesh.edge_length(m_inter.e)-m_inter.inter, 0);
                            else 
                                m_start_coord = Vec2(m_inter.inter, 0);
                        }
                        if (!m_min_iv->get_tau())
                            m_cur_e = m_min_e;
                        else
                            m_cur_e = m_min_e->Sym();
                        m_to_face = m_cur_e->Lface();
                        m_s = m_mesh.get_s_coord(m_min_iv);
                        m_ray = Ray<Vec2>(m_start_coord, m_s);

                        m_state = Test_s_state;
                    }
                }
                else
                    m_state = Tracing_ring_trace_state;
            }
        }
        break;
    }
    case Tracing_path_state:
    {
        Interval* next_iv = trace_ray_once(m_to_face, m_cur_e, m_ray, m_inter);
        if (next_iv != NULL) 
        {
            // store the path knot and coordinate
            store_path_knot(m_target_vertex, m_inter.e, m_inter.inter);

            // prepare for the next step tracing
            m_min_iv = next_iv;
            m_prev_e = m_min_e;
            if (m_min_iv->handle()->Org() == m_inter.e->Org())
            {
                if (!m_min_iv->get_tau())
                    m_start_coord = Vec2(m_inter.inter, 0);
                else 
                    m_start_coord = Vec2(m_mesh.edge_length(m_inter.e)-m_inter.inter, 0);
            } else {
                if (!m_min_iv->get_tau())
                    m_start_coord = Vec2(m_mesh.edge_length(m_inter.e)-m_inter.inter, 0);
                else 
                    m_start_coord = Vec2(m_inter.inter, 0);
            }
        }

        m_min_e = m_min_iv->handle();
        if (!m_min_iv->get_tau()) {
            m_cur_e = m_min_e;
        } else {
            m_cur_e = m_min_e->Sym();
        }
        m_to_face = m_cur_e->Lface();
        m_s = m_mesh.get_s_coord(m_min_iv);
        m_ray = Ray<Vec2>(m_start_coord, m_s);

        m_state = Test_s_state;
    }
    case Test_s_state:
    {
        if (FEQ(norm(m_start_coord-m_s), 0.0))
            m_state = Test_sigma_state;
        else 
            m_state = Tracing_path_state;
        break;
    }
    case Test_sigma_state:
    {
        m_next_vertex = get_interval_end_vertex(m_min_iv, m_start_coord);
        if (FEQ(m_min_iv->get_sigma(), 0.0))
            m_state = Complete_state;
        else 
            m_state = Tracing_ring_build_state;
        break;
    }
    case Complete_state: 
    {
        m_vindex++;
        m_state = Init_state;
    }
    }
}

void PathTracer::store_path_knot(int source_v, int knot_index, Handle eh)
{
    // store the encoded knot 
    PathKnotPair kp;
    HandleSet a_edges;
    Handle e = eh;
    if (e == NULL)
    {
        m_mesh.m_graph->collect_vertex_adj_edge(knot_index, a_edges);
        e = *a_edges.begin();
    }
    if (e == NULL) return;
    if (e->Org() == knot_index) 
        kp = make_pair(e, 0);
    else 
        kp = make_pair(e, m_mesh.edge_length(e));

    m_path_knot[source_v].push_back(kp);
}

void PathTracer::store_path_knot(int source_v, Handle e, double interval)
{
    // store the encoded knot
    PathKnotPair kp = make_pair(e, interval);
    m_path_knot[source_v].push_back(kp);
}

int PathTracer::get_next_vertex(Interval* iv)
{
    int vertex_id = -1;
    Handle e = iv->handle();
    if (FEQ(iv->get_d0(), 0.0))
        vertex_id = e->Org();
    else if (FEQ(iv->get_d1(), 0.0))
        vertex_id = e->Dest();

    return vertex_id;
}

void PathTracer::build_current_ring(int vid, Handle e, map<Handle, bool>& current_ring)
{
    HandleSet a_edges;
    m_mesh.m_graph->collect_vertex_adj_edge(vid, a_edges);
    
    flag_current_ring(current_ring, e);

    HandleSet::iterator hit = a_edges.begin();
    for ( ; hit != a_edges.end(); hit++)
    {
        if ((*hit) != e)
            current_ring[*hit] = false;
    }
}

Interval* PathTracer::get_one_interval(map<Handle, bool>& current_ring, Vec2& start_coord)
{
    Interval *min_iv = NULL;
    double min_dist(DBL_MAX);

    map<Handle, bool>::iterator hit = current_ring.begin();
    for ( ; hit != current_ring.end(); hit++) 
    {
        if ((*hit).second) continue;

        Handle e = (*hit).first;
        bool inverted = false;
        if (e->Org() > e->Dest()) {
            inverted = true;
            e = e->Sym();
        }
        const EdgeStruct& es = m_mesh.m_edge_map[e];
        if (!inverted) {
            double dist = es.get_distance_end0();
            if (dist < min_dist) {
                min_dist = dist;
                min_iv = es.front();
                if (!min_iv->get_tau())
                    start_coord = Vec2(0,0);
                else 
                    start_coord = Vec2(m_mesh.edge_length(e),0);
            }
        } else {
            double dist = es.get_distance_end1();
            if (dist < min_dist) {
                min_dist = dist;
                min_iv = es.back();
                if (!min_iv->get_tau())
                    start_coord = Vec2(m_mesh.edge_length(e),0);
                else 
                    start_coord = Vec2(0,0);
            }
        }
    }
    return min_iv;
}

void PathTracer::flag_current_ring(map<Handle, bool>& current_ring, Handle e)
{
    if (e != NULL) 
    {
        if (current_ring.find(e) != current_ring.end())
            current_ring[e] = true;
        if (current_ring.find(e->Sym()) != current_ring.end())
            current_ring[e->Sym()] = true;
    }
}

Interval* PathTracer::trace_ray_once(int to_face, Handle cur_e, Ray<Vec2>& ray, InterStruct& inter)
{
    id_t other_v_idx = m_mesh.other_vertex(to_face, cur_e->Org(), cur_e->Dest());
    Vec2 other_v = m_mesh.project_p(cur_e, to_face, m_mesh.m_vertex[other_v_idx]);
    double e_length = m_mesh.edge_length(cur_e);
    Interval* next_iv = NULL;

    // build the vindex coordinate map
    map<int, Vec2> coord_map;
    coord_map[cur_e->Org()]  = Vec2(0,0);
    coord_map[other_v_idx]   = other_v;
    coord_map[cur_e->Dest()] = Vec2(e_length, 0);

    bool intersected = m_mesh.intersect_face_border_loose(to_face, cur_e->Sym(), coord_map, ray, inter);

    if (intersected) 
        next_iv = get_interval_at(inter.e, inter.inter);

    return next_iv;
}

Interval* PathTracer::get_interval_at(Handle& eh, double inter)
{
    Handle e = eh;
    double interval = inter;
    if (e->Org() > e->Dest()) {
        e = e->Sym();
        interval = m_mesh.edge_length(e) - interval;
    }
    EdgeStruct& es = m_mesh.m_edge_map[e];
    return es.get_interval(interval);
}

int PathTracer::get_interval_end_vertex(Interval* iv, Vec2& start_coord) 
{
    Handle e = iv->handle();
    bool end0 = FEQ(iv->get_b0(), 0.0);
    bool end1 = FEQ(iv->get_b1(), m_mesh.edge_length(e));
    if (end0 && end1) {
        if (FEQ(start_coord[0], 0.0))
            if (!iv->get_tau())
                return e->Org();
            else 
                return e->Dest();
        else 
            if (!iv->get_tau())
                return e->Dest();
            else
                return e->Org();
    } else if (end0)
        return e->Org();

    return e->Dest();
}