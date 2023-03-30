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
    cout << "GeoTriMesh::compute_geodesic: " << selected_v << endl;
    reset_distance();

    if (is_first_propagation)
        init_propagation(selected_v);
    cout << "GeoTriMesh::compute_geodesic: " << selected_v << " " << m_vertex.size() << endl;

    // propagate intervals and compute distance field
    while (propagate_once());
    // trace the geodesic paths
    backtracing();
    cout << "GeoTriMesh::compute_geodesic - after backtracking: " << selected_v << " " << m_vertex.size() << endl;

    // output the results
    distance.resize(m_vertex.size());
    copy(geo_distance.begin(), geo_distance.end(), distance.begin());

    path.resize(m_vertex.size());
    copy(geo_path.begin(), geo_path.end(), path.begin());
}

bool GeoTriMesh::init_propagation(int selected_v)
{
    cout << "GeoTriMesh::init_propagation: " << selected_v << endl;
    Vec3& source_v = m_vertex[selected_v];

    // initialize the first elements in the priority queue
    HandleSet f_edges, a_edges;
    graph->collect_vertex_front_edge(selected_v, f_edges);
    graph->collect_vertex_adj_edge_with_sym(selected_v, a_edges);

    cout << "GeoTriMesh::init_propagation 2: " << selected_v << endl;
    HandleSet::iterator hit;
    
    // initialize intervals for the adjacent edges and insert them into edge map
    for (hit = a_edges.begin(); hit != a_edges.end(); hit++)
    {
        Interval* iv = create_init_interval(*hit, source_v);
        edge_map[(Handle) iv->handle()].push_back(iv);
    }

    cout << "GeoTriMesh::init_propagation 3: " << selected_v << endl;

    // initialize interval for the front edges, 
    // inser them into edge map, and insert them into priority queue
    for (hit = f_edges.begin(); hit != f_edges.end(); hit++)
    {
        Interval* iv = create_init_interval(*hit, source_v);
        edge_map[(Handle) iv->handle()].push_back(iv);
        i_queue.insert(iv);
    }

    cout << "GeoTriMesh::init_propagation 4: " << selected_v << endl;

    geo_distance[m_selected_vertex] = 0.f;
    for (EdgeStructMap::iterator eit = edge_map.begin(); eit != edge_map.end(); eit++)
        update_sfield((*eit).first, (*eit).second);

    is_first_propagation = false;

    cout << "GeoTriMesh::init_propagation 5: " << selected_v << endl;

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
        Handle to_e = (Handle) propag.altered.front()->handle();

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
    cout << "GeoTriMesh::update_sfield 1: " << endl;
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
    cout << "GeoTriMesh::update_sfield 2: " << endl;
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
    double elength = edge_length((Handle) iv->handle());
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
    double elength = edge_length((Handle) iv->handle());
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
    Handle e = (Handle) iv->handle();
    double dist = dvalue_at(iv, inter);

    return iv->sigma + dist;
}

double GeoTriMesh::dvalue_at(const Interval* iv, const double inter) const
{
    Handle e = (Handle) iv->handle();
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

void GeoTriMesh::merge_interval(Interval* new_iv, EdgeStruct& old_ivs, Propagation& propag)
{
	// merge the intersected intervals with iv
	Interval* v1 = new_iv;

	for (EdgeStruct::iterator oit = old_ivs.begin(); oit != old_ivs.end(); oit++)
	{
		// rename the intervals for clearity
		Interval* v0 = (*oit);

		if (!v0->does_intersect(v1))
			propag.insert_altered(v0); 
		else 
		{
			Range range = get_intersection_range(v0, v1);
			int compare_code = absolute_compare(v0, v1, range);

			switch (compare_code)
			{
			case -1:
			{
				// v0 has shorter distance all over the range
				if (DEBUG)
					cout << "case -1: v0 has shorter distance all over the range\n";
				absolute_merging_v0_shorter_v1(v0, v1, propag);
				break;
			}
			case 1:
			{
				// v1 has shorter distance all over the range
				if (DEBUG)
					cout << "case 1: v1 has shorter distance all over the range\n";
				absolute_merging_v1_shorter_v0(v0, v1, propag, range);
				break;
			}
			default:
			{
				// "sort" the two intervals
				if (DEBUG)
					cout << "case partial merging\n";
				partial_merging(v0, v1, propag, range);
			}
			}
		}
	}
	// fill gap (useful when v1 covers regions not covered by the existing intervals)
	// and clean really small intervals
	post_processing(propag, v1);
	propag.insert_deleted(v1);
}

Range GeoTriMesh::get_intersection_range(const Interval* v0, const Interval* v1) const
{
	Range rg;
	rg.first = FLT(v0->get_b0(), v1->get_b0()) ? v1->get_b0() : v0->get_b0();
	rg.second = FGT(v0->get_b1(), v1->get_b1()) ? v1->get_b1() : v0->get_b1();

	return rg;
}

short GeoTriMesh::absolute_compare(const Interval* v0, const Interval* v1, const Range& range) const
{
	double v0_0 = distance_at(v0, range.first);
	double v0_1 = distance_at(v0, range.second);
	double v1_0 = distance_at(v1, range.first);
	double v1_1 = distance_at(v1, range.second);

	if (!FGT(v0_0, v1_0) && !FGT(v0_1, v1_1))		return -1;	// v0 has shorter distance all over the range
	else if (!FLT(v0_0, v1_0) && !FLT(v0_1, v1_1))	return 1;	// v1 has shorter distance all over the range

	return 0;
}

void GeoTriMesh::absolute_merging_v0_shorter_v1(Interval* v0, Interval* v1, Propagation& propag)
{
	propag.insert_altered(v0);
}

void GeoTriMesh::absolute_merging_v1_shorter_v0(Interval* v0, Interval* v1, Propagation& propag, Range& rg)
{
	if (FLT(v0->get_b0(), rg.first))
		propag.insert_altered(v0->copy(v0->get_b0(), rg.first, v0->get_d0(), dvalue_at(v0, rg.first)));

	propag.insert_altered(v1->copy(rg.first, rg.second, dvalue_at(v1, rg.first), dvalue_at(v1, rg.second)));

	if (FGT(v0->get_b1(), rg.second))
		propag.insert_altered(v0->copy(rg.second, v0->get_b1(), dvalue_at(v0, rg.second), v0->get_d1()));

	propag.insert_deleted(v0);
}

double GeoTriMesh::get_px_value(Interval* v0, Interval* v1, Range& rg)
{
	// find the common point and update the intervals
	Vec2 s0 = get_s_coord_straight(v0);
	Vec2 s1 = get_s_coord_straight(v1);

	double alpha = s1[0] - s0[0];
	double beta = v1->get_sigma() - v0->get_sigma();
	double gamma = norm2(s0) - norm2(s1) - beta*beta;

	double A = alpha * alpha - beta * beta;
	double B = gamma * alpha + 2.f * s1[0] * beta * beta;
	double C = .25f * gamma * gamma - norm2(s1) * beta * beta;

	// it is stated in the paper that
	// it is a quadratic with a single solution within the range
	double px, delta = B*B - 4.0*A*C;
	if (delta <= 0)
		px = (-B) / (2.0*A);
	else if (delta > 0)
	{
		delta = sqrt(delta);
		px = (-B+delta) / (2.0*A);
		if (!GMath::within_range(px, rg))
			px = (-B-delta) / (2.0*A);
	}
	if (!GMath::within_range(px, rg)) {
		cout << "px: " << px << "   rg: " << rg.first << " " << rg.second << endl; 
	}
	assert(GMath::within_range(px, rg));

	return px;
}

void GeoTriMesh::partial_merging(Interval* v0, Interval* v1, Propagation& propag, Range& rg)
{
	double px = get_px_value(v0, v1, rg);

	unsigned char cond = 0;
	cond += test100(v0, v1, rg);
	cond += test010(v0, v1, rg);
	cond += test001(v0, v1, rg);

	switch (cond)
	{
	// for cases 0 and 1: v0->get_b0() == range.first
	//					  v0->get_b1() == range.second
	case 0:
		propag.insert_altered(v0->copy(v0->get_b0(), rg.first, v0->get_d0(), dvalue_at(v0, rg.first)));
		if (FLT(rg.first, px))
			propag.insert_altered(v1->copy(rg.first, px, dvalue_at(v1, rg.first), dvalue_at(v1, px)));
		if (FLT(px, v0->get_b1()))
			propag.insert_altered(v0->copy(px, v0->get_b1(), dvalue_at(v0, px), v0->get_d1()));
		break;
	case 1:
		if (FLT(v1->get_b0(), px))
			propag.insert_altered(v0->copy(v0->get_b0(), px, v0->get_d0(), dvalue_at(v0, px)));
		if (FLT(px, rg.second))
			propag.insert_altered(v1->copy(px, rg.second, dvalue_at(v1, px), dvalue_at(v1, rg.second)));
		break;
	// for cases 2 and 3: v1->get_b0() == range.first
	//					  v1->get_b1() == range.second
	case 2:
		propag.insert_altered(v0->copy(v0->get_b0(), rg.first, v0->get_d0(), dvalue_at(v0, rg.first)));
		if (FLT(rg.first, px))
			propag.insert_altered(v1->copy(rg.first, px, dvalue_at(v1, rg.first), dvalue_at(v1, px)));
		if (FLT(px, v0->get_b1()))
			propag.insert_altered(v0->copy(px, v0->get_b1(), dvalue_at(v0, px), v0->get_d1()));
		break;
	case 3:
		if (FLT(v0->get_b0(), px))
			propag.insert_altered(v0->copy(v0->get_b0(), px, v0->get_d0(), dvalue_at(v0, px)));
		if (FLT(px, rg.second))
			propag.insert_altered(v1->copy(px, rg.second, dvalue_at(v1, px), dvalue_at(v1, rg.second)));
		propag.insert_altered(v0->copy(rg.second, v0->get_b1(), dvalue_at(v0, rg.second), v0->get_d1()));
		break;
	// for cases 4 and 5: v0->get_b0() == range.first
	//					  v0->get_b1() == range.second
	case 4:
		if (FLT(rg.first, px))
			propag.insert_altered(v1->copy(rg.first, px, dvalue_at(v1, rg.first), dvalue_at(v1, px)));
		if (FLT(px, rg.second))
			propag.insert_altered(v0->copy(px, rg.second, dvalue_at(v0, px), dvalue_at(v0, rg.second)));
		break;
	case 5:
		if (FLT(rg.first, px))
			propag.insert_altered(v0->copy(rg.first, px, dvalue_at(v0, rg.first), dvalue_at(v0, px)));
		if (FLT(px, rg.second))
			propag.insert_altered(v1->copy(px, rg.second, dvalue_at(v1, px), dvalue_at(v1, rg.second)));
		break;
	// for cases 6 and 7: v0->get_b0() == range.first
	//					  v1->get_b1() == range.second
	case 6:
		if (FLT(rg.first, px))
			propag.insert_altered(v1->copy(rg.first, px, dvalue_at(v1, rg.first), dvalue_at(v1, px)));
		if (FLT(px, v0->get_b1()))
			propag.insert_altered(v0->copy(px, v0->get_b1(), dvalue_at(v0, px), v0->get_d1()));
		break;
	case 7:
		if (FLT(v0->get_b0(), px))
			propag.insert_altered(v0->copy(v0->get_b0(), px, v0->get_d0(), dvalue_at(v0, px)));
		if (FLT(px, rg.second))
			propag.insert_altered(v1->copy(px, rg.second, dvalue_at(v1, px), dvalue_at(v1, rg.second)));
		propag.insert_altered(v0->copy(rg.second, v0->get_b1(), dvalue_at(v0, rg.second), v0->get_d1()));
		break;
	}

	propag.insert_deleted(v0);
}

unsigned char GeoTriMesh::test100(Interval* v0, Interval* v1, Range& range)
{
	unsigned char res = (!FGT(range.first, v0->get_b0())) ? 4 : 0;
	return res;
}

unsigned char GeoTriMesh::test010(Interval* v0, Interval* v1, Range& range)
{
	unsigned char res = (FLT(range.second, v0->get_b1())) ? 2 : 0;
	return res;
}

unsigned char GeoTriMesh::test001(Interval* v0, Interval* v1, Range& range)
{
	unsigned char res = (FLT(distance_at(v0, range.first), distance_at(v1, range.first))) ? 1 : 0;
	return res;
}

void GeoTriMesh::post_processing(Propagation& propag, Interval* b)
{
	IntervalList& l = propag.altered;

	// fill the gaps
	if (!l.empty())
	{
		// copy the interval list to a temp list
		IntervalList temp_a(l.size());
		copy(l.begin(), l.end(), temp_a.begin());
		propag.clear_altered();

		double b_begin = b->get_b0();
		double b_end = b->get_b1();
		double a_end = temp_a.back()->get_b1();

		for (IntervalList::iterator ait = temp_a.begin(); ait != temp_a.end(); ait++)
		{
			Interval* a = *ait;
			double a_begin = a->get_b0();
			
			if (FGT(a_begin, b_begin))
				propag.insert_altered(b->copy(b_begin, a_begin, dvalue_at(b, b_begin), dvalue_at(b, a_begin)));

			l.push_back(a);
			b_begin = a->get_b1();
		}

		if (FGT(b_end, a_end))
			propag.insert_altered(b->copy(a_end, b_end, dvalue_at(b, a_end), dvalue_at(b, b_end)));
	}

	// clean up very small intervals
	int lsize = l.size();
	if (lsize > 1) 
	{
		IntervalList& d = propag.deleted;
		vector<bool> remains(lsize, true);
		IntervalList temp_l(lsize);
		copy(l.begin(), l.end(), temp_l.begin());

		Interval* cur = temp_l.front();
		IntervalList::iterator iit;
		int i=0;
		for (iit = temp_l.begin()+1; iit != temp_l.end(); iit++, i++)
		{
			Interval* next = (*iit);

			if (cur->is_very_small())
			{
				next->set_range(cur->get_b0(), next->get_b1(), dvalue_at(next, cur->get_b0()), next->get_d1());
				remains[i] = false;
			}
			cur = next;
		}

		l.clear();
		for (iit = temp_l.begin(), i=0; iit != temp_l.end(); iit++, i++)
		{
			if (remains[i])	l.push_back(*iit);
			else			propag.insert_deleted(*iit);
		}
	}
}

void GeoTriMesh::propagate_interval(Interval* iv, EdgeStructMap& new_ivs)
{
	Handle cur_e, from_e = (Handle) iv->handle();
	double from_e_length = edge_length(from_e);
	id_t from_face, to_face;
	Vec2 b0_pos, b1_pos;
	bool tau = iv->get_tau();
	
	if (graph->is_infinite_face(from_e->Sym()->Lface()) || graph->is_infinite_face(from_e->Lface()))
		return;
	
	if (!tau) {
		// the interval is not inverted
		cur_e = from_e;
		from_face = from_e->Lface();
		to_face = from_e->Sym()->Lface();
		b0_pos = Vec2(iv->get_b0(), 0.0);
		b1_pos = Vec2(iv->get_b1(), 0.0);
	} else {
		// the interval is inverted
		cur_e = from_e->Sym();
		from_face = from_e->Sym()->Lface();
		to_face = from_e->Lface();
		b0_pos = Vec2(from_e_length-iv->get_b1(), 0.0);
		b1_pos = Vec2(from_e_length-iv->get_b0(), 0.0);
	}

	Vec2 s = get_s_coord(iv);
	last_step.i0_pt = Vec3();
	last_step.i1_pt = Vec3();
	last_step.ray0 = Ray<Vec2>();
	last_step.ray1 = Ray<Vec2>();

	if (FEQ(s[1], 0.0)) {
		if ( s[0] <= 0.0 ) {
			double cur_sigma = iv->get_sigma() - s[0];
			cover_diamond_shape(iv, cur_e, m_vertex[cur_e->Org()], cur_sigma, from_e_length, new_ivs, true);
		}
		else if ( s[0] >= from_e_length ) {
			double cur_sigma = iv->get_sigma() + s[0] - from_e_length;
			cover_diamond_shape(iv, cur_e, m_vertex[cur_e->Dest()], cur_sigma, from_e_length, new_ivs, false);
		}
	}
	else 
	{
		int other_v_idx = other_vertex(to_face, from_e->Org(), from_e->Dest());
		Vec2 other_v = project_p(cur_e, to_face, m_vertex[other_v_idx]);
		last_step.other_v = unproject_p(cur_e, to_face, other_v);

		InterStruct b0_inter, b1_inter;
		Ray<Vec2> b0_ray(s, b0_pos);
		Ray<Vec2> b1_ray(s, b1_pos);

		// build the vindex coordinate map
		map<int, Vec2> coord_map;
		coord_map[cur_e->Org()]  = Vec2(0,0);
		coord_map[other_v_idx]    = other_v;
		coord_map[cur_e->Dest()] = Vec2(from_e_length, 0);

		// compute the intersections
		bool b0_intersected = intersect_face_border(to_face, cur_e, coord_map, b0_ray, b0_inter);
		bool b1_intersected = intersect_face_border(to_face, cur_e, coord_map, b1_ray, b1_inter);

		if (!b0_intersected || !b1_intersected) {
			last_step.i0_pt = (b0_inter.e == 0) ? Vec3(0) : edge_point(b0_inter.e, b0_inter.inter);
			last_step.i1_pt = (b1_inter.e == 0) ? Vec3(0) : edge_point(b1_inter.e, b1_inter.inter);
		} else {
			// case 2
			if (b0_inter.e != b1_inter.e)
				propagate_case2(b0_inter, b1_inter, norm(s-other_v), iv, new_ivs);

			// case 3 - additional on the left
			else if (is_iv_b0_end(tau, iv, from_e_length) && b0_inter.e == cur_e->Dnext())
				propagate_case3_left(b0_inter, b1_inter, cur_e, s, iv, new_ivs);
									  
			// case 3 - additional on the right
			else if (is_iv_b1_end(tau, iv, from_e_length) && b0_inter.e == cur_e->Oprev())
				propagate_case3_right(b0_inter, b1_inter, cur_e, s, iv, new_ivs);
			
			// case 1 - most simple propagation
			else 
				propagate_case1(b0_inter, b1_inter, iv, new_ivs);		

			last_step.i0_pt = edge_point(b0_inter.e, b0_inter.inter);
			last_step.i1_pt = edge_point(b1_inter.e, b1_inter.inter);
		}
		last_step.ray0 = b0_ray;
		last_step.ray1 = b1_ray;
	}

	last_step.cur_iv = iv;
	last_step.from_face = from_face;
	last_step.to_face = to_face;
	last_step.cur_e = cur_e;
	last_step.s = s;
}

void GeoTriMesh::propagate_case1(const InterStruct& b0_it, const InterStruct& b1_it, Interval* iv, 
								 EdgeStructMap& new_ivs)
{
	// the propagated interval
	if (FGT(b1_it.inter, b0_it.inter))
	{
		Interval* prop = create_interval(b0_it.e, b0_it.inter, b1_it.inter, b0_it.t, b1_it.t, iv->get_sigma());
		new_ivs[(Handle) prop->handle()].insert_interval(prop);
	}
}

void GeoTriMesh::propagate_case2(const InterStruct& b0_it, const InterStruct& b1_it, double mid_dist,
								 Interval* iv, EdgeStructMap& new_ivs)
{
	double sigma = iv->get_sigma();
	if (FGT(edge_length(b0_it.e), b0_it.inter))
	{
		Interval* i0 = create_interval(b0_it.e, b0_it.inter, edge_length(b0_it.e), b0_it.t, mid_dist, sigma);
		new_ivs[(Handle) i0->handle()].insert_interval(i0);
	}
	if (FGT(b1_it.inter, 0.f))
	{
		Interval* i1 = create_interval(b1_it.e, 0.f, b1_it.inter, mid_dist, b1_it.t, sigma);
		new_ivs[(Handle) i1->handle()].insert_interval(i1);
	}
}

void GeoTriMesh::propagate_case3_left(const InterStruct& b0_it, const InterStruct& b1_it, Handle from_e, 
									  const Vec2& s, Interval* iv, EdgeStructMap& new_ivs)
{
	Handle prop_e = b0_it.e;
	Handle gap_e;
	Vec2 pseudo_s;
	double from_e_length = edge_length(from_e);

	gap_e = from_e->Oprev();
//	pseudo_s = (iv->get_tau()) ? Vec2(from_e_length,0) : Vec2(0,0);
	pseudo_s = Vec2(0,0);

	double cur_sigma = iv->get_sigma() + norm(s-pseudo_s);
	double gap_e_length = edge_length(gap_e);

	// the propagated interval
	if (FGT(b1_it.inter, b0_it.inter))
	{
		Interval* prop = create_interval(prop_e, b0_it.inter, b1_it.inter, b0_it.t, b1_it.t, iv->get_sigma());
		new_ivs[(Handle) prop->handle()].insert_interval(prop);
	}
	// create two more intervals with pseudo source
	if (FGT(b0_it.inter, 0.f))
	{
		Interval* partial = create_interval(prop_e, 0.f, b0_it.inter, gap_e_length, b0_it.t-norm(s-pseudo_s), cur_sigma);
		new_ivs[(Handle) partial->handle()].insert_interval(partial);
	}
	Interval* full = create_interval(gap_e, 0.f, gap_e_length, 0, gap_e_length, cur_sigma);
	new_ivs[(Handle) full->handle()].insert_interval(full);
}

void GeoTriMesh::propagate_case3_right(const InterStruct& b0_it, const InterStruct& b1_it, Handle from_e, 
									   const Vec2& s, Interval* iv, EdgeStructMap& new_ivs)
{
	Handle prop_e = b0_it.e;
	Handle gap_e = from_e->Dnext();
	Vec2 pseudo_s;
	double from_e_length = edge_length(from_e);

//	pseudo_s = (iv->get_tau()) ? Vec2(0,0) : Vec2(from_e_length,0);
	pseudo_s = Vec2(from_e_length,0);

	double gap_e_length = edge_length(gap_e);
	double prop_e_length = edge_length(b1_it.e);
	double cur_sigma = iv->sigma + norm(s-pseudo_s);

	// the propagated interval
	if (FGT(b1_it.inter, b0_it.inter))
	{
		Interval* prop = create_interval(prop_e, b0_it.inter, b1_it.inter, b0_it.t, b1_it.t, iv->get_sigma());
		new_ivs[(Handle) prop->handle()].insert_interval(prop);
	}
	// create two more gap intervals with pseudo source
	if (FGT(prop_e_length, b1_it.inter))
	{
		Interval* partial = create_interval(prop_e, b1_it.inter, prop_e_length, b1_it.t-norm(s-pseudo_s), gap_e_length, cur_sigma);
		new_ivs[(Handle) partial->handle()].insert_interval(partial);
	}
	Interval* full = create_interval(gap_e, 0.f, gap_e_length, gap_e_length, 0, cur_sigma);
	new_ivs[(Handle) full->handle()].insert_interval(full);
}

void GeoTriMesh::cover_diamond_shape(Interval* iv, Handle cur_e, Vec3& source, double cur_sigma, double e_length, 
									 EdgeStructMap& new_ivs, bool start)
{
	Handle aux_e;
	if (!graph->is_infinite_face(cur_e->Lface())) 
	{
		aux_e = cur_e->Lnext();
		fill_diamond_gap(aux_e, cur_sigma, source, new_ivs);
		aux_e = cur_e->Lprev();	
		fill_diamond_gap(aux_e, cur_sigma, source, new_ivs);
	}
	if (!graph->is_infinite_face(cur_e->Sym()->Lface())) 
	{
		aux_e = cur_e->Oprev();
		fill_diamond_gap(aux_e, cur_sigma, source, new_ivs);
		aux_e = cur_e->Dnext();
		fill_diamond_gap(aux_e, cur_sigma, source, new_ivs);
	}
}

void GeoTriMesh::fill_diamond_gap(Handle cur_e, double cur_sigma, Vec3& source, EdgeStructMap& new_ivs)
{
	bool inverted = (cur_e->Org() > cur_e->Dest());
	Handle e = !inverted ? cur_e : cur_e->Sym();

	EdgeStruct& es = edge_map[e];
	if (es.empty()) 
	{
		new_ivs[e].insert_interval(create_interval(cur_e, cur_sigma, source));
	}
	else 
	{
		double b0_ext = es.get_smallest_b0(), b1_ext = es.get_largest_b1();
		double e_length = edge_length(e);
		// fill the small extrem
		if (FGT(b0_ext, 0))
			new_ivs[e].insert_interval(create_interval(e, cur_sigma, 0, b0_ext, source, inverted));
		// fill the intermediate gaps
		int esize = es.size();
		for (int i=1; i<esize; i++)
		{
			Interval* prev = es[i-1];
			Interval* next = es[i];
			if (FEQ(prev->get_b1(), next->get_b0(), 1e-4))
				new_ivs[e].insert_interval(create_interval(e, cur_sigma, prev->get_b1(), next->get_b0(), source, inverted));
		}
		// fill the big extrem
		if (FLT(b1_ext, e_length))
			new_ivs[e].insert_interval(create_interval(e, cur_sigma, b1_ext, e_length, source, inverted));
	}
}

bool GeoTriMesh::intersect_face_border(int face, Handle from_e, map<int,Vec2>& coord_map, Ray<Vec2> ray, InterStruct& inter)
{
	Handle adj_e = from_e->Sym();
	Handle e = adj_e->Lnext(), min_e;
	double min_t = FLT_MAX;
	bool intersected = false;
	Vec2 p0, p1, point, min_point;

	for (int i=0; i<2; i++)
	{
		// test the intersection for the other edges of the face
		p0 = coord_map[e->Org()];
		p1 = coord_map[e->Dest()];

		double t = ray.segment_intersect(p0, p1, point);
		if (ray.has_infinite_intersection())
		{
			Vec2& pt = (norm2(ray.p-p0) < norm2(ray.p-p1)) ? p0 : p1;
			double temp_t = norm(pt-ray.p);
			if (temp_t < min_t) 
			{
				min_t = norm(pt-ray.p);
				min_e = e;
				min_point = pt;
				intersected = true;
			}
		}
		else if (ray.has_one_intersection() && t > 0.0 && t < min_t)
		{
			min_t = t;
			min_e = e;
			min_point = point;
			intersected = true;
		}
		e = e->Lnext();
	}
//	assert(intersected);

	if (!intersected) {
		Handle e = adj_e->Lnext();
		if (DEBUG) cout << "adj e: " << *e << endl;
		int test_v = e->Org();
		if (close_enough(coord_map[test_v], ray))
		{
			inter.t = norm(coord_map[test_v]-ray.p);
			inter.e = e;
			inter.inter = 0.0;
			intersected = true;
		}
		e = e->Lnext();
		if (DEBUG) cout << "adj e: " << *e << endl;
		test_v = e->Dest();
		if (!intersected && close_enough(coord_map[test_v], ray))
		{
			inter.t = norm(coord_map[test_v]-ray.p);
			inter.e = e;
			inter.inter = edge_length(e);
			intersected = true;
		}
	}
	else
	{
		p0 = coord_map[min_e->Org()];
		p1 = coord_map[min_e->Dest()];
		inter.e = min_e;
		inter.inter = edge_length(min_e)*edge_interval_unit(p0, p1, min_point);
		inter.t = min_t;
	} 
	return intersected;
}


bool GeoTriMesh::intersect_face_border_loose(int face, Handle from_e, map<int,Vec2>& coord_map, Ray<Vec2> ray, InterStruct& inter)
{
	Handle adj_e = from_e->Sym();
	Handle e = adj_e->Lnext(), min_e;
	double min_t = FLT_MAX;
	bool intersected = false;
	Vec2 p0, p1, point, min_point;

	for (int i=0; i<2; i++)
	{
		// test the intersection for the other edges of the face
		p0 = coord_map[e->Org()];
		p1 = coord_map[e->Dest()];

		double t = ray.segment_intersect(p0, p1, point);
		if (ray.has_infinite_intersection())
		{
			Vec2& pt = (norm2(ray.p-p0) < norm2(ray.p-p1)) ? p0 : p1;
			double temp_t = norm(pt-ray.p);
			if (temp_t < min_t) 
			{
				min_t = norm(pt-ray.p);
				min_e = e;
				min_point = pt;
				intersected = true;
			}
		}
		else if (ray.has_one_intersection() && t>0.0 && t < min_t)
		{
			min_t = t;
			min_e = e;
			min_point = point;
			intersected = true;
		}
		e = e->Lnext();
	}

	if (!intersected) {
		Handle e = adj_e->Lnext();
		int test_v = e->Org();
		if (close_enough(coord_map[test_v], ray.p))
		{
			inter.t = 0.0;
			inter.e = e;
			inter.inter = 0.0;
			intersected = true;
		}
		e = e->Lnext();
		test_v = e->Dest();
		if (!intersected && close_enough(coord_map[test_v], ray.p))
		{
			inter.t = 0.0;
			inter.e = e;
			inter.inter = edge_length(e);
			intersected = true;
		}
	}
	else
	{
		p0 = coord_map[min_e->Org()];
		p1 = coord_map[min_e->Dest()];
		inter.e = min_e;
		inter.inter = edge_length(min_e)*edge_interval_unit(p0, p1, min_point);
		inter.t = min_t;
	} 
	return intersected;
}

bool GeoTriMesh::close_enough(const Vec2& p0, const Vec2& p1)
{
	if (DEBUG) double temp = norm2(p0-p1);
	return FEQ(norm2(p0-p1), 0.0);
}

bool GeoTriMesh::close_enough(const Vec2& p, const Ray<Vec2>& ray)
{
	double a = ray.d[1], b = -ray.d[0], c = - (a*ray.p[0] + b*ray.p[1]);
	double dist = a*p[0] + b*p[1] + c;
	return FEQ(dist, 0.0);
}