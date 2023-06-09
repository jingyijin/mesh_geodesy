/************************************************************************

  Common code for ball-based rotation controllers.

  Code originally part of LibQSlim; 
  Placed in this directory to keep the library files free of OpenGL dependencies.

  $Id: baseball.cxx 427 2004-09-27 04:45:31Z garland $

 ************************************************************************/

#include "baseball.hpp"

#include <sstream>

namespace gfx
{

Baseball::Baseball()
{
    curquat = Quat::ident();

    trans=0.0;
    ctr=0.0;
    radius=1;
}

void Baseball::apply_transform()
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glTranslated(trans[0], trans[1], trans[2]);
    glTranslated(ctr[0], ctr[1], ctr[2]);

    const Mat4 M=unit_quat_to_matrix(curquat);
    glMultMatrixd(M);

    glTranslated(-ctr[0], -ctr[1], -ctr[2]);
}

void Baseball::unapply_transform()
{
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void Baseball::write(std::ostream& out)
{
    out << "baseball ";
    out << curquat << " " << trans << " " << ctr << " " << radius << std::endl;
}

void Baseball::read(std::istream& in)
{
    std::string name;

    in >> name;
    in >> curquat >> trans >> ctr >> radius;
}

} // namespace gfx
