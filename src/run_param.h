#pragma once
#include <stdbool.h>

typedef double REAL;

#define NMESH_X (128)
#define NMESH_V (128)

#define NGREEN (NMESH_X/2+1)

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

#define NOUTPUT_MAX (128)

struct run_param {
  int nstep;
  REAL tnow, tend;
  REAL dtime;

  REAL xmax, xmin;
  REAL vmax, vmin;

  REAL delta_x, delta_v;

  int output_indx;
  REAL output_timing[NOUTPUT_MAX];

  char modelname[128];

  bool grav;  // true for gravitational interaction 
};
