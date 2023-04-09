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
#include "general_math.hpp"
#include "heap.hpp"

enum TraversalState { Discovered, Finished };

/**
 * @class Interval
 * @brief Represents a narrow band interval associated with a mesh edge.
 * The narrow band interval is used in the propagation of distances in the triangle mesh.
 * The class stores the mesh edge (halfedge) associated with the interval, its endpoints,
 * scalar values at the end points, the pseudo-source distance, side of the (pseudo) source
 * and the accumulated error in approximation. The class also defines functions to access
 * and modify its attributes, as well as to compute the minimum distance, to check for
 * intersection and to merge with another interval.
 */
class Interval : public MxHeapable
{
public:
    typedef Halfedge<id_t>::Handle Handle; /**< Handle to a halfedge */

public:
    Handle e;       /**< Handle to the halfedge associated with the interval */
    double b0, b1;  /**< Interval endpoints */
    double d0, d1;  /**< Scalar values at the endpoints */
    double sigma;   /**< Pseudo-source distance */
    bool tau;       /**< Side of the (pseudo) source:
                     * positive if the source is on e->LFace()
                     * negative if the source is on e->Sym()->LFace()
                     */
    double xi;      /**< Accumulated error in approximation */

    TraversalState state;   /**< State of the interval:
                             * Discovered if the interval is in the narrow band
                             * Finished if the interval is not in the narrow band
                             */

public:
    /**
     * @brief Constructor for an interval with given endpoint values and source attributes.
     * @param e_ Handle to the Halfedge this interval is associated with.
     * @param b0_ Left endpoint value of the interval.
     * @param b1_ Right endpoint value of the interval.
     * @param d0_ Scalar value at the left endpoint.
     * @param d1_ Scalar value at the right endpoint.
     * @param sigma_ Pseudo-source distance.
     * @param tau_ Side of the (pseudo) source. Positive if the source is on e->LFace(), negative if the source is on e->Sym()->LFace().
     */
    Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, bool tau_)
        : MxHeapable()
    {
        e = e_;
        b0 = b0_; b1 = b1_;
        d0 = d0_; d1 = d1_;
        sigma = sigma_;
        tau = tau_;
        xi = 0.0;
        state = Discovered;
        compute_min_dist();
    }

    /**
     * @brief Constructor for an interval with given endpoint values, source attributes, and traversal state.
     * @param e_ Handle to the Halfedge this interval is associated with.
     * @param b0_ Left endpoint value of the interval.
     * @param b1_ Right endpoint value of the interval.
     * @param d0_ Scalar value at the left endpoint.
     * @param d1_ Scalar value at the right endpoint.
     * @param sigma_ Pseudo-source distance.
     * @param tau_ Side of the (pseudo) source. Positive if the source is on e->LFace(), negative if the source is on e->Sym()->LFace().
     * @param state_ Traversal state of the interval.
     */    
    Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, bool tau_, TraversalState state_)
        : MxHeapable()
    {
        e = e_;
        b0 = b0_; b1 = b1_;
        d0 = d0_; d1 = d1_;
        sigma = sigma_;
        tau = tau_;
        xi = 0.0;
        state = state_;
        compute_min_dist();
    }
    /**
     * @brief Constructor for an interval with given endpoint values, source attributes, traversal state, and accumulated error in approximation.
     * @param e_ Handle to the Halfedge this interval is associated with.
     * @param b0_ Left endpoint value of the interval.
     * @param b1_ Right endpoint value of the interval.
     * @param d0_ Scalar value at the left endpoint.
     * @param d1_ Scalar value at the right endpoint.
     * @param sigma_ Pseudo-source distance.
     * @param xi_ Accumulated error in approximation.
     * @param tau_ Side of the (pseudo) source. Positive if the source is on e->LFace(), negative if the source is on e->Sym()->LFace().
     * @param state_ Traversal state of the interval.
     */    
    Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, double xi_, bool tau_, TraversalState state_)
        : MxHeapable()
    {
        e = e_;
        b0 = b0_; b1 = b1_;
        d0 = d0_; d1 = d1_;
        sigma = sigma_;
        tau = tau_;
        xi = xi_;
        state = state_;
        compute_min_dist();
    }

    inline Handle handle() const { return e; }          /**< Get the handle of the interval */
    inline double key() const { return heap_key(); }    /**< Get the heap key of the interval */
    inline double get_b0() const { return b0; }         /**< Get the interval endpoint b0 */
    inline double get_b1() const { return b1; }         /**< Get the interval endpoint b1 */
    inline double get_d0() const { return d0; }         /**< Get the scalar value at the interval endpoint b0 */
    inline double get_d1() const { return d1; }         /**< Get the scalar value at the interval endpoint b1 */
    inline double get_sigma() const { return sigma; }   /**< Get the pseudo-source distance */
    inline bool get_tau() const { return tau; }         /**< Get the side of the (pseudo) source */
    inline double get_xi() const { return xi; }         /**< Get the accumulated error in approximation */
    inline TraversalState get_state() const { return state; } /**< Get the traversal state of the interval */

    inline void set_b0(double b0_) { b0 = b0_; } /**< Set the value of b0 attribute */
    inline void set_b1(double b1_) { b1 = b1_; } /**< Set the value of b1 attribute */
    inline void set_d0(double d0_) { d0 = d0_; } /**< Set the value of d0 attribute */
    inline void set_d1(double d1_) { d1 = d1_; } /**< Set the value of d1 attribute */
    inline void set_xi(double xi_) { xi = xi_; } /**< Set the value of xi attribute */
    inline void set_sigma(double sigma_) { sigma = sigma_; } /**< Set the value of sigma attribute */

    inline bool is_discovered() const { return state == Discovered; } /**< Check if interval is discovered. */
    inline bool is_finished() const { return state == Finished; }     /**< Check if interval is finished. */
    inline void set_discovered() { state = Discovered; }              /**< Set interval as discovered. */
    inline void set_finished() { state = Finished; }                  /**< Set interval as finished. */

    /**
     * @brief Computes the minimum distance value from the pseudo-source to one of the end points.
     * @return Minimum distance value.
     */
    inline double min_dist() const
    {
        return sigma + ((d0 < d1) ? d0 : d1);
    }

    /**
     * @brief Computes the minimum distance value and sets it as the heap key for the Interval instance.
     */
    inline void compute_min_dist()
    {
        heap_key(-min_dist());
    }
    inline double get_d0_dist() const { return sigma + get_d0(); } /**< Returns the distance to the pseudo-source at d0. */
    inline double get_d1_dist() const { return sigma + get_d1(); } /**< Returns the distance to the pseudo-source at d1. */

    inline bool is_b0_end(double e_value) const { return FEQ(b0, e_value); } /**< Checks if b0 is equal to the given value e_value. */
    inline bool is_b1_end(double e_value) const { return FEQ(b1, e_value); } /**< Checks if b1 is equal to the given value e_value. */
    /**
     * @brief Checks if this interval intersects with the given Interval instance.
     * @param iv Interval instance to check for intersection with this interval.
     * @return True if the two intervals intersect, false otherwise.
     */
    inline bool does_intersect(Interval* iv) const
    {
        return ( FLT(get_b0(), iv->get_b1()) && FGT(get_b1(), iv->get_b0()) );
    }
    /**
     * @brief Checks if this interval covers the whole segment between two given end points.
     * @param end0 The first end point.
     * @param end1 The second end point.
     * @return True if this interval covers the whole segment between the two given end points, false otherwise.
     */
    inline bool cover_whole(double end0, double end1) const
    {
        return ( FEQ(get_b0(), end0) && FEQ(get_b1(), end1) );
    }
    /**
     * @brief Checks if the interval is very small, i.e., if the difference between its two endpoints is close to zero.
     * @return True if the interval is very small, false otherwise.
     */
    inline bool is_very_small() const 
    {
        return FEQ(get_b0()-get_b1(), 0.0, 1e-4);
    }

    /**
     * @brief Creates a copy of this Interval instance with new values for its attributes.
     * @param new_b0 New value for the b0 attribute of the new Interval instance.
     * @param new_b1 New value for the b1 attribute of the new Interval instance.
     * @param new_d0 New value for the d0 attribute of the new Interval instance.
     * @param new_d1 New value for the d1 attribute of the new Interval instance.
     * @return A pointer to the new Interval instance.
     */
    inline Interval* copy(double new_b0, double new_b1, double new_d0, double new_d1) 
    {
        return new Interval(e, new_b0, new_b1, new_d0, new_d1, sigma, tau, state);
    }
    /**
     * @brief Updates the values of the attributes of this Interval instance with the new given values.
     * @param new_b0 New value for the b0 attribute of this Interval instance.
     * @param new_b1 New value for the b1 attribute of this Interval instance.
     * @param new_d0 New value for the d0 attribute of this Interval instance.
     * @param new_d1 New value for the d1 attribute of this Interval instance.
     */
    inline void set_range(double new_b0, double new_b1, double new_d0, double new_d1)
    {
        b0 = new_b0; b1 = new_b1;
        d0 = new_d0; d1 = new_d1;
    }
    /**
     * Merges two Intervals by updating the attributes of this Interval based on the attributes of Interval a and a new point new_s
     * @param a The Interval to merge with this Interval
     * @param new_s A new point to update the Interval's attributes based on
     * @param new_sigma The new value of sigma for the Interval
     * @param new_xi The new value of xi for the Interval
     */
    inline void merge(Interval* a, Vec2& new_s, double new_sigma, double new_xi)
    {
        this->set_b0(a->get_b0());
        this->set_d0(norm(new_s-Vec2(a->get_b0(), 0)));
        this->set_d1(norm(new_s-Vec2(this->get_b1(), 0)));
        this->set_sigma(new_sigma);
        this->set_xi(new_xi);
    }

    /**
     * Overloads the output stream operator to print the Interval's attributes
     * @param out The output stream to print to
     * @param iv The Interval to print
     * @return A reference to the output stream
     */
    friend std::ostream& operator<<(std::ostream& out, const Interval& iv)
    { 
        return out << "e: " << (*iv.handle()) << " b0: " << iv.get_b0() << " b1: " << iv.get_b1() 
                << " d0: " << iv.get_d0() << " d1: " << iv.get_d1()
                << " sigma: " << iv.get_sigma() << " tau: " << iv.get_tau() << " finished: " << iv.is_finished(); 
    }
};


#endif