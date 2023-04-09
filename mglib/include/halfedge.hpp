#ifndef HALF_EDGE_H_INCLUDED 
#define HALF_EDGE_H_INCLUDED

/************************************************************************
 * File description: Half-edge data structure for mesh generation
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/

#include <assert.h>
#include <ostream>
#include <vector>

typedef unsigned int iid_t;
typedef std::pair<iid_t, iid_t> idpair_t;
typedef std::vector<iid_t> idlist_t;

/**
 * @brief A simple structure to represent a 2D polygon.
 */
struct Polygon : public idlist_t
{
    Polygon(size_t N) : idlist_t(N) {}
    Polygon(const idlist_t& p) : idlist_t(p) {}

    int dim() const { return 2; }
};

/**
 * @brief A class representing half-edge data structure.
 * @tparam _Vertex The vertex type.
 * @tparam _Face The face type (default: iid_t).
 */
template<class _Vertex, class _Face=iid_t>
class Halfedge
{
public:
    typedef Halfedge *Handle;
    typedef _Vertex Vertex;
    typedef _Face Face;
    enum { NoFace = -1 };

public:
    Handle lnext, lprev, sym;
    Vertex org;
    Face lface;

public:
    /**
     * @brief Construct a half-edge with the given vertex and face.
     * @param v The vertex of the half-edge.
     * @param f The face of the half-edge (default: NoFace).
     */
    Halfedge(Vertex v, Face f=NoFace)
        : lnext(nullptr), lprev(nullptr), sym(nullptr), org(v), lface(f) {}
    /**
     * @brief Destructor for the half-edge.
     */
    ~Halfedge() { lnext=lprev=sym=nullptr; }

    Handle Lnext() const { return lnext; } 
    Handle Lprev() const { return lprev; }
    Handle Sym() const { return sym; }

    Handle Onext() const { return Lprev()->Sym(); }
    Handle Oprev() const { return Sym()->Lnext(); }
    Handle Dnext() const { return Sym()->Lprev(); }
    Handle Dprev() const { return Lnext()->Sym(); }
    Handle Rnext() const { return Oprev()->Sym(); }
    Handle Rprev() const { return Dnext()->Sym(); }

    Vertex& Org()		 { return org; }
    Vertex  Org() const  { return org; }

    // we use Lnext here, rather than Sym because Lnext will
    // always be defined, whereas Sym may not be.
    Vertex& Dest()		 { return Lnext()->Org(); }
    Vertex  Dest() const { return Lnext()->Org(); }

    // Halfedges are oriented from e[0] -> e[1]
    Vertex& operator[](int i)		{ return i ? Dest() : Org(); }
    Vertex  operator[](int i) const { return i ? Dest() : Org(); }

    // Each edge maintains an ID for the face on its left
    Face& Lface()		{ return lface; }
    Face  Lface() const { return lface; }

public:
    /**
     * @brief Create a face from a given polygon.
     * @tparam Polygon The polygon type.
     * @param P The input polygon.
     * @param f The face identifier (default: NoFace).
     * @return The handle of the first half-edge of the face.
     */
    template<class Polygon>
    static Handle create_face(Polygon& P, Face f=NoFace)
    {
        int i;
        const int N = P.size();
        std::vector<Handle> edges(N);

        for (i=0; i<N; i++) edges[i] = new Halfedge(P[i], f);
        for (i=0; i<N-1; i++) edges[i]->lnext = edges[i+1];
        for (i=1; i<N; i++) edges[i]->lprev = edges[i-1];

        edges[N-1]->lnext = edges[0];
        edges[0]->lprev = edges[N-1];

        return edges[0];
    }

    /**
     * @brief Connect two half-edges with a symmetric relationship.
     * @param e1 The first half-edge.
     * @param e2 The second half-edge.
     */
    static void paste(Handle e1, Handle e2)
    {
        assert(e1->Sym() == nullptr);
        assert(e2->Sym() == nullptr);
        assert(e1->Org() == e2->Dest());

        e1->sym = e2;
        e2->sym = e1;
    }
    
    /**
     * @brief Disconnect the symmetric relationship between two half-edges.
     * @param e The half-edge to be cut.
     */
    static void cut(Handle e)
    {
        if (e->sym)
        {
            e->sym->sym = nullptr;
            e->sym = nullptr;
        }
    }

   /**
     * @brief Stream out the half edge.
     * @param e The half-edge to be streamed.
     */
    template<class T>
    friend inline std::ostream& operator<<(std::ostream& out, const Halfedge<T>& e)
    { 
        return out << e.Org() << " " << e.Dest(); 
    }
};


#endif