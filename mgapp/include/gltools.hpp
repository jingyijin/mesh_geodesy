#ifndef GFXGLTOOLS_INCLUDED // -*- C++ -*-
#define GFXGLTOOLS_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  Handy functions for common OpenGL tasks
  Code originally part of LibQSlim; 
  it is pulled out of gfx library to keep gfx clean of OpenGL dependencies.

  $Id: gltools.cxx 446 2005-06-18 13:58:15Z garland $

 ************************************************************************/

#include "vec4.hpp"

#include <cassert>

namespace gfx
{

using std::cerr;
using std::endl;

GLuint opengl_pick_nil = (~0);
GLuint opengl_pick_zmax = (~0);

void begin_opengl_pick(int *where, double radius, GLuint *buffer, int size)
{
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);

    glSelectBuffer(size, buffer);
    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(opengl_pick_nil);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();             // Save the current transformation
    glLoadIdentity();
    gluPickMatrix(where[0], vp[3] - where[1], radius, radius, vp);
}

GLuint complete_opengl_pick(GLuint *buffer)
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();              // get rid of the pick matrix
    glFlush();

    GLint nhits = glRenderMode(GL_RENDER);
    GLuint *hit = NULL;
    GLuint hit_nnames = 0;
    GLuint zmin = opengl_pick_zmax;
    GLuint *ptr = buffer;

    for(int i=0; i<nhits; i++)
    {
        GLuint nnames   = *ptr++;
        GLuint cur_zmin = *ptr++;
        /* GLuint cur_zmax = */ *ptr++;

        if( cur_zmin < zmin )
        {
            zmin = cur_zmin;
            hit = ptr;
	    hit_nnames = nnames;
        }
        ptr+=nnames;
    }


    buffer[0] = hit_nnames;
    if( hit )
    {
	for(int k=0; k<hit_nnames; k++)
	    buffer[k+1] = hit[k];

	return *hit;
    }
    else
	return opengl_pick_nil;
}

void camera_lookat(const Vec3& min, const Vec3& max, double aspect)
{
    Vec3 up(0, 1, 0);
    double fovy = 60.0;

    Vec3 at = (max + min)/2.0;         // look at the center of the bbox
    double radius = norm(max - at);    // radius of a bounding sphere
    double d = 3*radius / tan(fovy * M_PI/180.0);

    Vec3 from = at;
    from[2] += d;

    double znear = d/20;
    double zfar = 10*d;

    glMatrixMode(GL_PROJECTION);
    gluPerspective(fovy, aspect, znear, zfar);


    glMatrixMode(GL_MODELVIEW);
    gluLookAt(from[0], from[1], from[2],
	      at[0], at[1], at[2],
	      up[0], up[1], up[2]);
}


int unproject_pixel(int *pixel, double *world, double z)
{
    GLdouble modelMatrix[16];
    GLdouble projMatrix[16];
    GLint viewport[4];

    glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
    glGetIntegerv(GL_VIEWPORT, viewport);

    // Notice that we have to correct the y pixel coordinate.  GL
    // assigns the origin to the lower left corner, while FLTK assigns
    // the origin to the upper left corner.
    return gluUnProject(pixel[0], viewport[3]-pixel[1], z,
			modelMatrix, projMatrix, viewport,
			world, world+1, world+2);
}

} // namespace gfx

// GFXGLTOOLS_INCLUDED
#endif