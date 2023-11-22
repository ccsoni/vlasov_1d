#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "run_param.h"
#include "vlasov_1d.h"
#include "prototype.h"

#define DF(ix,iv) df[(iv)+NMESH_V*(ix)]

void setup_output_timing(struct run_param *this_run, REAL dt_output)
{
  this_run->output_indx = 0;

  for(int32_t i=0;i<NOUTPUT_MAX;i++) {
    this_run->output_timing[i] = dt_output*(REAL)(i);
  }
}

int main(int argc, char **argv)
{
  struct run_param this_run;

  REAL *df, *dens, *pot;

  df   = (REAL *)malloc(sizeof(REAL)*NMESH_X*NMESH_V);
  dens = (REAL *)malloc(sizeof(REAL)*NMESH_X);
  pot  = (REAL *)malloc(sizeof(REAL)*NMESH_X);

  REAL dtime;

  setup_IC(df, &this_run);
  setup_output_timing(&this_run, 0.1);

  calc_dens_field(df, dens, &this_run);
  calc_poten(dens, pot, &this_run);
  dtime = calc_dtime(&this_run, pot);

  output_diag(df, &this_run, dens, pot);

  while(this_run.tnow < this_run.tend) {
    // advance along x by a half time step
    integrate_x(df, 0.5*dtime, &this_run);

    if(this_run.grav) {
      // advance along v by a full time step
      calc_dens_field(df, dens, &this_run);
      calc_poten(dens, pot, &this_run);
      integrate_v(df, pot, dtime, &this_run);
    }

    // advance along x by a half time step    
    integrate_x(df, 0.5*dtime, &this_run);
    
    this_run.tnow += dtime;
    this_run.nstep++;
    output_diag(df, &this_run, dens, pot);
    output_df_in_run(df, &this_run);

    dtime = calc_dtime(&this_run, pot);
  }

  free(dens);
  free(pot);
  free(df);
}

