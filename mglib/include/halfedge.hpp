#ifndef HALF_EDGE_H_INCLUDED 
#define HALF_EDGE_H_INCLUDED

#include "cell_types.hpp"
#include <assert.h>
#include <ostream>

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
    Halfedge(Vertex v, Face f=NoFace)
        : lnext(nullptr), lprev(nullptr), sym(nullptr), org(v), lface(f) {}

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

    static void paste(Handle e1, Handle e2)
    {
        assert(e1->Sym() == nullptr);
        assert(e2->Sym() == nullptr);
        assert(e1->Org() == e2->Dest());

        e1->sym = e2;
        e2->sym = e1;
    }
    
    static void cut(Handle e)
    {
        if (e->sym)
        {
            e->sym->sym = nullptr;
            e->sym = nullptr;
        }
    }
};

template<class T>
inline std::ostream& operator<<(std::ostream& out, const Halfedge<T>& e)
{ return out << e.Org() << " " << e.Dest(); }

#endif