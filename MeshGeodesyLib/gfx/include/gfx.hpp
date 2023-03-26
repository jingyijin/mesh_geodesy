#ifndef GFX_INCLUDED // -*- C++ -*-
#define GFX_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  Main header file for the libgfx graphics library.

  $Id: gfx.h 455 2005-08-17 18:10:25Z garland $

 ************************************************************************/

#include <cstdlib>
#include <cmath>
#include <climits>
#include <iostream>


#ifndef M_PI
#  define M_PI 3.14159265358979323846264338327950288419716939937510582097494459
#endif


////////////////////////////////////////////////////////////////////////
//
namespace gfx
{

#if defined(HAVE_RANDOM)
  inline double random1() { return (double)random() / (double)LONG_MAX; }
  inline char   random_byte() { return (char)(random() & 0xff); }
#else
  inline double random1() { return (double)rand() / (double)RAND_MAX; }
  inline char   random_byte() { return (char)(rand() & 0xff); }
#endif

const double FEQ_EPS = 1e-6;
const double FEQ_EPS2 = 1e-12;

inline bool  FEQ(double a, double b, double e=FEQ_EPS)  {return fabs(a-b)<e;}
inline bool FEQ2(double a, double b, double e=FEQ_EPS2) {return fabs(a-b)<e;}


////////////////////////////////////////////////////////////////////////
//
//
//

#define TIMING(t, cmd) { t=get_cpu_time(); cmd; t=get_cpu_time() - t; }
extern double get_cpu_time();

} // namespace gfx

#ifndef GFX_NAMESPACE
using namespace gfx;
#endif

// GFX_INCLUDED
#endif
