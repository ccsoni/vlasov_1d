#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>

#include "run_param.h"
#include "vlasov_1d.h"
#include "prototype.h"

#define DF(ix,iv) df[(iv)+NMESH_V*(ix)]

void setup_IC_rect(REAL *df, struct run_param *this_run)
{
  this_run->xmin = 0.0;
  this_run->xmax = 1.0;
  this_run->delta_x = (this_run->xmax-this_run->xmin)/(REAL)NMESH_X;

  this_run->vmin = -1.0;
  this_run->vmax = 1.0;
  this_run->delta_v = (this_run->vmax-this_run->vmin)/(REAL)NMESH_V;

  for(int32_t ix=0;ix<NMESH_X;ix++) {
    for(int32_t iv=0;iv<NMESH_V;iv++) {
      REAL x = this_run->xmin + this_run->delta_x*((REAL)ix+0.5);
      REAL v = this_run->vmin + this_run->delta_v*((REAL)iv+0.5);

      if(0.25 < x && x<0.75 && -0.5<v && v<0.5) {
	DF(ix,iv) = 1.0;
      }else{
	DF(ix,iv) = 0.0;
      }
    }
  }

  this_run->tend = 5.0;
  this_run->grav = false;
  
  sprintf(this_run->modelname, "%s", "rect");  
}

void setup_IC_two_stream(REAL *df, struct run_param *this_run)
{
  int32_t wave_num = 1;
  REAL delta = 0.01;

  this_run->xmin = 0.0;
  this_run->xmax = 1.0;
  this_run->delta_x = (this_run->xmax-this_run->xmin)/(REAL)NMESH_X;

  this_run->vmin = -2.0;
  this_run->vmax = 2.0;
  this_run->delta_v = (this_run->vmax-this_run->vmin)/(REAL)NMESH_V;

  REAL kwave = 2.0*M_PI/(this_run->xmax-this_run->xmin)*(REAL)wave_num;
  REAL kjeans = kwave/0.4;
  REAL beta = 1.0;

  REAL vel_disp = 4.0*M_PI/SQR(kjeans);
  REAL v_t = sqrt(vel_disp)*5.0;

  for(int32_t ix=0;ix<NMESH_X;ix++) {
    REAL x = this_run->xmin + this_run->delta_x*((REAL)ix+0.5);
    for(int32_t iv=0;iv<NMESH_V;iv++) {
      REAL v  = this_run->vmin + this_run->delta_v*((REAL)iv+0.5);
      REAL _erf_p =
	erf((v+0.5*this_run->delta_v-0.5*v_t)/sqrt(2.0*vel_disp))
        -erf((v-0.5*this_run->delta_v-0.5*v_t)/sqrt(2.0*vel_disp));
      REAL _erf_m =
        erf((v+0.5*this_run->delta_v+0.5*v_t)/sqrt(2.0*vel_disp)/beta)
        -erf((v-0.5*this_run->delta_v+0.5*v_t)/sqrt(2.0*vel_disp)/beta);

      DF(ix,iv) =
        0.5*(_erf_p + _erf_m)*(1.0+delta*cos(kwave*x))/this_run->delta_v;
    }
  }

  this_run->tend = 5.0;
  this_run->grav = true;

  sprintf(this_run->modelname, "%s", "two_stream");
}

void setup_IC_collapse(REAL *df, struct run_param *this_run)
{
  int32_t wave_num = 2;
  REAL delta = 0.01; // density fluctuation

  this_run->xmin = 0.0;
  this_run->xmax = 1.0;
  this_run->delta_x = (this_run->xmax-this_run->xmin)/(REAL)NMESH_X;

  this_run->vmin = -1.0;
  this_run->vmax = 1.0;
  this_run->delta_v = (this_run->vmax-this_run->vmin)/(REAL)NMESH_V;

  REAL kwave = 2.0*M_PI/(this_run->xmax-this_run->xmin)*(REAL)wave_num;
  REAL kjeans = kwave/0.5;  
  REAL vel_disp = 4.0*M_PI/SQR(kjeans);

  for(int32_t ix=0;ix<NMESH_X;ix++) {
    REAL x = this_run->xmin + this_run->delta_x*((REAL)ix+0.5);
    for(int32_t iv=0;iv<NMESH_V;iv++) {
      REAL v  = this_run->vmin + this_run->delta_v*((REAL)iv+0.5);
      REAL _erf = erf((v+0.5*this_run->delta_v)/sqrt(2.0*vel_disp))
        -erf((v-0.5*this_run->delta_v)/sqrt(2.0*vel_disp));
      DF(ix,iv) = _erf*0.5*(1.0+delta*cos(kwave*x))/this_run->delta_v;
    }
  }

  this_run->tend = 5.0;
  this_run->grav = true;

  sprintf(this_run->modelname, "%s", "collapse");
}

void setup_IC(REAL *df, struct run_param *this_run)
{
  this_run->tnow = 0.0;
  this_run->nstep = 0;
  
  //setup_IC_rect(df, this_run);
  setup_IC_collapse(df, this_run);
  //setup_IC_two_stream(df, this_run);
}
