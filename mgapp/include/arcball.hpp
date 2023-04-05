#ifndef GFXARCBALL_INCLUDED // -*- C++ -*-
#define GFXARCBALL_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  Arcball rotation control.
  Code originally part of LibQSlim; 
  it is pulled out of gfx library to keep gfx clean of OpenGL dependencies.
  
  $Id: arcball.h 427 2004-09-27 04:45:31Z garland $

 ************************************************************************/

#include "baseball.hpp"

namespace gfx
{

/**
@brief Arcball class, a sub-class of Baseball for arcball interaction.
*/

class Arcball : public Baseball
{
private:
    Vec2 ball_ctr;             // Center of arcball in pixels
    double ball_radius;        // Radius of arcball in pixels

    Quat q_now, q_down, q_drag;	// Quaternions describing rotation
    Vec3 v_from, v_to;	        // Mouse drag vectors

    bool is_dragging;          // True if mouse is being dragged

protected:
    /**
    * @brief Project a point on the arcball sphere from a 2D screen point.
    * @param p Point to be projected
    * @return The point on the sphere that the 2D point is projected onto
    */
    Vec3 proj_to_sphere(const Vec2&);
    /**
     * @brief Update the state of the arcball based on the current quaternion
     */
    void update();


public:
    /**
    * @brief Default constructor for Arcball class
    */
    Arcball();

    /**
     * @brief Update the animation of the Arcball
     */
    virtual void update_animation();
    /**
     * @brief Handle the mouse down event for the Arcball
     * @param where The 2D location of the mouse down event
     * @param which Which mouse button is being used for this event
     * @return True if the mouse down event is handled successfully
     */    
    virtual bool mouse_down(int *where, int which);
    /**
     * @brief Handle the mouse up event for the Arcball
     * @param where The 2D location of the mouse up event
     * @param which Which mouse button is being used for this event
     * @return True if the mouse up event is handled successfully
     */
    virtual bool mouse_up(int *where, int which);
    /**
     * @brief Handle the mouse drag event for the Arcball
     * @param where The current 2D location of the mouse drag event
     * @param last The previous 2D location of the mouse drag event
     * @param which Which mouse button is being used for this event
     * @return True if the mouse drag event is handled successfully
     */    
    virtual bool mouse_drag(int *where, int *last, int which);

    /**
     * @brief Apply the current transformation of the Arcball
     */
    virtual void apply_transform();
    /**
     * @brief Get the current transformation of the Arcball
     * @param c Center of the arcball
     * @param t Translation of the arcball
     * @param q Rotation of the arcball
     */    
    virtual void get_transform(Vec3 & c, Vec3 &t, Quat & q);
    /**
     * @brief Set the current transformation of the Arcball
     * @param c Center of the arcball
     * @param t Translation of the arcball
     * @param q Rotation of the arcball
     */    
    virtual void set_transform(const Vec3 & c, const Vec3 & t, const Quat & q); 

    /**
     * @brief Write the Arcball transformation to an output stream
     * @param os The output stream to write to
     */
    virtual void write(std::ostream&);
    /**
     * @brief Read the Arcball transformation from an input stream
     * @param is The input stream to read from
     */
    virtual void read(std::istream&);
};

} // namespace gfx

// GFXARCBALL_INCLUDED
#endif
