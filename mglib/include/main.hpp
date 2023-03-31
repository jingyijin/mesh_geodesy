#ifndef MAIN_INCLUDED
#define MAIN_INCLUDED

#include "gfx.hpp"
#include "vec3.hpp"
#include <glog/logging.h>

#include <vector>

using namespace std;

const double FEQ_EPS_LOOSE = 1e-3;

inline bool FLT(double a, double b, double e=FEQ_EPS)  { return (b-a)>e;}
inline bool FGT(double a, double b, double e=FEQ_EPS)  { return (a-b)>e;}
inline bool FLT_LOOSE(double a, double b, double e=FEQ_EPS_LOOSE)  { return (b-a)>e;}
inline bool FGT_LOOSE(double a, double b, double e=FEQ_EPS_LOOSE)  { return (a-b)>e;}

extern double current_time, prev_time;

#endif