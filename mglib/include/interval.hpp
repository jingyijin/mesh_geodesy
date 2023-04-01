#ifndef INTERVAL_H_INCLUDED
#define INTERVAL_H_INCLUDED

/************************************************************************
 * File description: Interval structure for the wave front algorithm
 *
 * Author: Jingyi Jin
 * Version: 0.1
 * Date: 4/1/2023
 ************************************************************************/



#include "halfedge.hpp"
#include "main.hpp"
#include "heap.hpp"

enum TraversalState { Discovered, Finished };

class Interval : public MxHeapable
{
public:
    typedef Halfedge<id_t>::Handle Handle;

public:
    Handle e;
    double b0, b1;	// interval endpoints
    double d0, d1;	// scalar values at the end points
    double sigma;	// pseudo-source distance
    bool tau;		// side of the (pseudo) source
                    // positive if the source is on e->LFace()
                    // negative if the source is on e->Sym()->LFace()
    double xi;		// accumulated error in approximation

    TraversalState state;	// traversal state

public:
    Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, bool tau_);
    Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, bool tau_, TraversalState state_);
    Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, double xi_, bool tau_, TraversalState state_);

    // accessing attributes
    Handle handle() const { return e; }
    double key() const    { return heap_key(); }
    double get_b0() const { return b0; }
    double get_b1() const { return b1; }
    double get_d0() const { return d0; }
    double get_d1() const { return d1; }
    double get_sigma() const { return sigma; }
    bool get_tau() const { return tau; }
    double get_xi() const { return xi; }
    TraversalState get_state() const { return state; }

    // setting values
    void set_b0(double b0_) { b0 = b0_; }
    void set_b1(double b1_) { b1 = b1_; }
    void set_d0(double d0_) { d0 = d0_; }
    void set_d1(double d1_) { d1 = d1_; }
    void set_xi(double xi_) { xi = xi_; }
    void set_sigma(double sigma_) { sigma = sigma_; }

    bool is_discovered() const { return state == Discovered; }
    bool is_finished() const { return state == Finished; }
    void set_discovered() { state = Discovered; }
    void set_finished() { state = Finished; }

    double min_dist() const;
    void compute_min_dist();
    double get_d0_dist() { return sigma + get_d0(); }
    double get_d1_dist() { return sigma + get_d1(); }

    bool is_b0_end(double e_value) const { return FEQ(b0, e_value); }
    bool is_b1_end(double e_value) const { return FEQ(b1, e_value); }
    
    bool does_intersect(Interval* iv) const;
    bool cover_whole(double end0, double end1) const;
    bool is_very_small() const;

    Interval* copy(double new_b0, double new_b1, double new_d0, double new_d1);
    void set_range(double new_b0, double new_b1, double new_d0, double new_d1);
    void merge(Interval* b, Vec2& new_s, double new_sigma, double new_xi);
};

inline std::ostream& operator<<(std::ostream& out, const Interval& iv)
{ 
    return out << "e: " << (*iv.handle()) << " b0: " << iv.get_b0() << " b1: " << iv.get_b1() 
               << " d0: " << iv.get_d0() << " d1: " << iv.get_d1()
               << " sigma: " << iv.get_sigma() << " tau: " << iv.get_tau() << " finished: " << iv.is_finished(); 
}

#endif