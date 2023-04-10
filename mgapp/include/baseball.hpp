#ifndef GFXBASEBALL_INCLUDED // -*- C++ -*-
#define GFXBASEBALL_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  Common base class for ball-based rotators (e.g., Trackball & Arcball).
  Code originally part of LibQSlim; 
  it is pulled out of gfx library to keep gfx clean of OpenGL dependencies.
  
  $Id: baseball.h 443 2005-06-14 00:53:40Z garland $

 ************************************************************************/

#include "quat.hpp"

#include <GL/gl.h>
#include <GL/glu.h>

namespace gfx
{

/**
@brief The Baseball class represents a controller for an object with a bounding sphere.
*/
class Baseball
{
public:
    Vec3 ctr;       // Describes bounding sphere of object
    double radius;  // Describes bounding sphere of object

    Quat curquat;   // Current rotation of object
    Vec3 trans;     // Current translation of object

public:
    /**
    * @brief Default constructor for the Baseball class.
    */
    Baseball();
    /**
     * @brief Virtual destructor for the Baseball class.
     */
    virtual ~Baseball() {}

    /**
     * @brief Required initialization method for setting the bounding sphere of the object.
     * @tparam T The type of vector elements (float, double, etc.).
     * @param v The center of the bounding sphere.
     * @param r The radius of the bounding sphere.
     */
    template<class T>
    void bounding_sphere(const TVec3<T>& v, T r) { ctr=v; radius=r; }

    /**
     * @brief Virtual method for updating the animation.
     */
    virtual void update_animation() = 0;
    /**
     * @brief Virtual method for handling mouse down events.
     * @param where Pointer to an array of integers representing the position of the mouse.
     * @param which Integer representing the type of mouse button clicked.
     * @return Boolean indicating if the event was handled.
     */
    virtual bool mouse_down(int *where, int which) = 0;
    /**
     * @brief Virtual method for handling mouse up events.
     * @param where Pointer to an array of integers representing the position of the mouse.
     * @param which Integer representing the type of mouse button clicked.
     * @return Boolean indicating if the event was handled.
     */    
    virtual bool mouse_up(int *where, int which) = 0;
    /**
     * @brief Virtual method for handling mouse drag events.
     * @param where Pointer to an array of integers representing the position of the mouse.
     * @param last Pointer to an array of integers representing the previous position of the mouse.
     * @param which Integer representing the type of mouse button clicked.
     * @return Boolean indicating if the event was handled.
     */
    virtual bool mouse_drag(int *where, int *last, int which) = 0;

    /**
     * @brief Interface method for applying the appropriate transformation during drawing.
     */
    virtual void apply_transform();
    /**
     * @brief Interface method for reversing the transformation applied during drawing.
     */    
    virtual void unapply_transform();

    /**
     * @brief Interface method for writing the transform to a stream.
     * @param[out] os The output stream.
     */
    virtual void write(std::ostream&);
    /**
     * @brief Interface method for reading the transform from a stream.
     * @param[in] is The input stream.
     */
    virtual void read(std::istream&);
};

} // namespace gfx

// GFXBASEBALL_INCLUDED
#endif
