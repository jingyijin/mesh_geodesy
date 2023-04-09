#ifndef MANIFOLD_H_INCLUDED
#define MANIFOLD_H_INCLUDED

/************************************************************************
 * File description: Manifold graph data structure for mesh generation
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include "halfedge.hpp"

#include <map>

/**
 * @brief A template class for representing a manifold graph with a specified topology.
 * @tparam _Topology The topology of the manifold graph.
 */
template<class _Topology>
class ManifoldGraph
{
public:
    typedef _Topology Topology;
    typedef typename Topology::Cell Cell;
    typedef typename Topology::CellList CellList;
    typedef typename Topology::VertexList VertexList;

    typedef id_t CellID;
    typedef id_t VertexID;
    typedef Halfedge<id_t> Edge;
    typedef typename Edge::Handle Handle;
    typedef std::set<Handle> HandleSet;
    typedef std::map<idpair_t, Edge*> EdgeMap;

public:
    CellList& m_cells;        /**< A reference to the list of cells in the manifold graph. */

    EdgeMap m_open_edges;     /**< A map of open edges in the manifold graph. */
    size_t m_halfedge_count;  /**< The number of halfedges in the manifold graph. */
    size_t m_vertex_count;    /**< The number of vertices in the manifold graph. */

    std::vector<Handle> m_infinite_cell;          /**< A vector of handles to the infinite cells (boundary) in the manifold graph. */
    std::set<id_t> m_infinite_cell_id;            /**< A set of IDs of the infinite cells in the manifold graph. */
    std::map<id_t, HandleSet> m_edges_of_vertices;/**< A map of edges of vertices in the manifold graph. */

public:
    /**
     * @brief Constructor for the ManifoldGraph class.
     *
     * @param nvert The number of vertices in the manifold graph.
     * @param C A reference to the list of cells in the manifold graph.
     */
    ManifoldGraph(size_t nvert, CellList& C)
        : m_cells(C), m_vertex_count(nvert), m_halfedge_count(0)
    {
        link_finite_elements();
        link_infinite_elements();
    }
    /**
     * @brief Destructor for the ManifoldGraph class.
     */
    ~ManifoldGraph()
    {
        // clear open edges
        typename EdgeMap::iterator e_it, e_end = m_open_edges.end();
        for(e_it = m_open_edges.begin(); e_it != e_end; e_it++)
            delete e_it->second;
        // clear edges of vertices
        typename std::map<id_t, HandleSet>::iterator v_it, v_end = m_edges_of_vertices.end();
        for(v_it = m_edges_of_vertices.begin(); v_it != v_end; v_it++)
            v_it->second.clear();
    }

    /**
     * @brief Link finite elements of the cell list.
     * It iterates over all the cells in the list, creates halfedges for each side, and links them to form a triangular mesh.
     */
    void link_finite_elements()
    {
        typename CellList::iterator e_it, e_end = m_cells.end();
        id_t v_ids[4], v0_id, v1_id, fid=0;
        Handle e_hs[3];
        for(e_it = m_cells.begin(); e_it != e_end; e_it++, fid++)
        {
            v_ids[0] = (*e_it)[0];	v_ids[1] = (*e_it)[1];
            v_ids[2] = (*e_it)[2];	v_ids[3] = (*e_it)[0];

            for(int i=0; i<3; i++)
            {
                v0_id = v_ids[i];	v1_id = v_ids[i+1];
                e_hs[i] = new Edge(v0_id, fid);
                m_edges_of_vertices[v0_id].insert(e_hs[i]);

                idpair_t p(v0_id, v1_id);
                idpair_t q(v1_id, v0_id);

                if( m_open_edges.find(q) == m_open_edges.end() )
                    m_open_edges[p] = e_hs[i];
                else
                {
                    Handle e2 = m_open_edges[q];
                    assert( e_hs[i]->Sym() == NULL );
                    assert( e2->Sym() == NULL );
                    assert( e_hs[i]->Org() == e2->Dest() );

                    e_hs[i]->sym = e2;
                    e2->sym = e_hs[i];

                    m_open_edges.erase(q);
                }
            }

            e_hs[0]->lnext = e_hs[1];	e_hs[1]->lnext = e_hs[2];
            e_hs[2]->lnext = e_hs[0];

            e_hs[1]->lprev = e_hs[0];	e_hs[2]->lprev = e_hs[1];
            e_hs[0]->lprev = e_hs[2];
        }
    }

    /**
     * @brief Links infinite elements in the mesh.
     * This function links infinite elements in the mesh by iterating over open edges,
     * creating new infinite elements and connecting them properly. It also stores
     * infinite elements in the infinite_cell vector and updates the infinite_cell_id set.
     */
    void link_infinite_elements(void)
    {
        //now, collect infinite elements
        EdgeMap::iterator o_e_it, o_e_end = m_open_edges.end();
        id_t infinite_e_id = -2;
        id_t v0_id, v1_id;
        unsigned int size, i;
        for (o_e_it = m_open_edges.begin(); o_e_it != o_e_end; o_e_it++)
        {
            Handle e_h = o_e_it->second;
            assert(e_h->Sym() == NULL);
            v0_id = e_h->Dest();	
            v1_id = e_h->Org();
            Handle inf_e_h = new Edge(v0_id, infinite_e_id);
            m_edges_of_vertices[v0_id].insert(inf_e_h);
            inf_e_h->sym = e_h;	e_h->sym = inf_e_h;
        }

        while ( (o_e_it = m_open_edges.begin()) != m_open_edges.end() )
        {
            std::vector<Handle> edges;
            Handle e_h = o_e_it->second;
            while (true)
            {
                v0_id = e_h->Dest();	v1_id = e_h->Org();
                Handle inf_e_h = e_h->Sym();
                size = m_open_edges.erase(idpair_t(v1_id, v0_id));
                if(size == 0)	break;
                edges.push_back(inf_e_h);

                e_h = NULL;
                HandleSet::iterator v_e_it, v_e_end;
                v_e_end = m_edges_of_vertices[v1_id].end();
                for (v_e_it = m_edges_of_vertices[v1_id].begin(); v_e_it != v_e_end; v_e_it++)
                {
                    Handle tmp_eh = *v_e_it;
                    if(tmp_eh->Lnext() == NULL)
                    {
                        e_h = tmp_eh->Sym();
                        break;
                    }
                }
            }

            size = edges.size() - 1;
            assert(size >= 2);
            for(i=0; i<size; i++)
            {
                edges[i]->lnext = edges[i+1];
                edges[i+1]->lprev = edges[i];
                edges[i]->lface = infinite_e_id;
            }

            edges[size]->lnext = edges[0];
            edges[0]->lprev = edges[size];
            edges[size]->lface = infinite_e_id;
            m_infinite_cell.push_back(edges[0]);
            m_infinite_cell_id.insert(infinite_e_id);
            infinite_e_id--;
        }
    }

public:
    // Emulation of selected CellGraph methods
    /**
     * @brief Collects adjacent edges of a vertex.
     * @param vid Vertex ID.
     * @param ering Set to store adjacent edges.
     */    
    void collect_vertex_adj_edge(const VertexID vid, HandleSet& ering)
    {
        assert(vertex_exist(vid));
        ering.clear();
        HandleSet& es = m_edges_of_vertices[vid];
        for (HandleSet::iterator eit = es.begin(); eit != es.end(); eit++) 
        {
            const Handle eh = *eit;
            if (is_infinite_face(eh->Lface())) continue;
            ering.insert(eh);
        }
    }

    /**
     * @brief Collects adjacent edges of a vertex along with their symmetric edges.
     * @param vid Vertex ID.
     * @param ering Set to store adjacent edges with symmetric edges.
     */
    void collect_vertex_adj_edge_with_sym(const VertexID vid, HandleSet& ering)
    {
        assert(vertex_exist(vid));
        ering.clear();
        HandleSet& es = m_edges_of_vertices[vid];
        for (HandleSet::iterator eit = es.begin(); eit != es.end(); eit++) 
        {
            const Handle eh = *eit;
            if (!is_infinite_face(eh->Lface()))
                ering.insert(eh);
            else 
                ering.insert(eh->Sym());
        }
    }

    /**
     * @brief Collects front edges of a vertex.
     *
     * @param vid Vertex ID.
     * @param ering Set to store front edges.
     */    
    void collect_vertex_front_edge(const VertexID vid, HandleSet& ering)
    {
        assert(vertex_exist(vid));
        ering.clear();
        HandleSet& es = m_edges_of_vertices[vid];
        for (HandleSet::iterator eit = es.begin(); eit != es.end(); eit++) 
        {
            const Handle eh = *eit;
            if (is_infinite_face(eh->Lface())) continue;
            ering.insert(eh->Lnext());
        }
    }

    /**
     * @brief Collects vertices forming a ring around a given vertex.
     *
     * @param vid Vertex ID.
     * @param vring Set to store the ring of vertices.
     */
    void collect_vertex_ring(const VertexID vid, std::set<VertexID>& vring)
    {
        assert(vertex_exist(vid));
        vring.clear();
        HandleSet& ering = m_edges_of_vertices[vid];
        for (HandleSet::iterator eit = ering.begin(); eit != ering.end(); eit++)
        {
            const Handle eh = *eit;
            if (is_infinite_face(eh->Lface())) continue;
            assert(eh->Org() == vid);
            vring.insert(eh->Dest());
        }
    }

    /**
     * @brief Collects faces forming a ring around a given vertex.
     *
     * @param vid Vertex ID.
     * @param fring Set to store the ring of faces.
     */
    void collect_face_ring(const VertexID vid, std::set<CellID>& fring)
    {
        assert(vertex_exist(vid));
        fring.clear();
        HandleSet& ering = m_edges_of_vertices[vid];
        for (HandleSet::iterator eit = ering.begin(); eit != ering.end(); eit++)
        {
            const Handle eh = *eit;
            if (is_infinite_face(eh->Lface())) continue;
            assert(eh->Org() == vid);
            fring.insert(eh->Lface());
        }
    }

    /**
     * @brief Checks if a vertex exists in the mesh.
     *
     * @param vid Vertex ID.
     * @return True if the vertex exists, false otherwise.
     */
    bool vertex_exist(const VertexID vid) const
    {	return vid >= 0 && vid < m_vertex_count;	}

    /**
     * @brief Checks if a face is an infinite face.
     *
     * @param cid Cell ID (face ID).
     * @return True if the face is infinite, false otherwise.
     */
    bool is_infinite_face(const CellID cid) const
    {	return cid < 0 || m_infinite_cell_id.find(cid) != m_infinite_cell_id.end();	 }
};

#endif