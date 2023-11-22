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

void output_diag(REAL *df, struct run_param *this_run, REAL *dens, REAL *pot)
{
  static FILE *diag_file;
  static int diag_file_open = 0;

  if(diag_file_open == 0) {
    //diag_file = fopen("vlasov_1d.diag","w");
    diag_file = stdout;
    diag_file_open = 1;
    fprintf(diag_file, "# istep    tnow    dtime    mass   KE    PE    TE \n");
  }

  REAL KE, PE;

  // recalculate the density field because it is over-written in calc_poten().
  calc_dens_field(df, dens, this_run);

  KE = PE = 0.0;
  for(int ix=0;ix<NMESH_X;ix++) {
    PE += 0.5*dens[ix]*pot[ix]*this_run->delta_x;
    for(int iv=0;iv<NMESH_V;iv++) {
      REAL vel = this_run->vmin + this_run->delta_v*((REAL)iv+0.5);
      REAL v_up = vel + 0.5*this_run->delta_v;
      REAL v_lo = vel - 0.5*this_run->delta_v;
      KE += 0.5*DF(ix,iv)*(CUBE(v_up)-CUBE(v_lo))/3.0;
    }
  }
  KE *= this_run->delta_x;

  REAL mass;

  mass = 0.0;
  for(int ix=0;ix<NMESH_X;ix++) {
    mass += dens[ix]*this_run->delta_x;
  }

  fprintf(diag_file, "%d %17.9e %14.6e %17.9e %17.9e %17.9e %17.9e\n",
          this_run->nstep, this_run->tnow, this_run->dtime, mass, KE, PE, KE+PE);

}

void output_df(REAL *df, struct run_param *this_run, char *filename)
{
  FILE *output_fp;

  output_fp = fopen(filename,"w");
#if 0
  fwrite(this_run, sizeof(struct run_param), 1, output_fp);
  fwrite(df, sizeof(REAL), NMESH_X*NMESH_V, output_fp);
#else
  fprintf(output_fp, "#%d %d \n", NMESH_X, NMESH_V);
  fprintf(output_fp, "#%12.4e %12.4e \n", this_run->xmin, this_run->xmax);
  fprintf(output_fp, "#%12.4e %12.4e \n", this_run->vmin, this_run->vmax);
  for(int32_t ix=0;ix<NMESH_X;ix++) {
    for(int32_t iv=0;iv<NMESH_V;iv++) {
      fprintf(output_fp, "%d %d %14.6e\n",ix, iv, DF(ix,iv));
    }
    fprintf(output_fp, "\n");
  }
#endif

  fclose(output_fp);
}

void input_df(REAL *df, struct run_param *this_run, char *filename)
{
  FILE *input_fp;

  input_fp = fopen(filename, "r");

#if 0
  fread(this_run, sizeof(struct run_param), 1, input_fp);
  fread(df, sizeof(REAL), NMESH_X*NMESH_V, input_fp);
#else
  int32_t io_err;
  int32_t nmesh_x, nmesh_y;
  io_err = fscanf(input_fp, "%d %d", &nmesh_x, &nmesh_y);
  assert(nmesh_x == NMESH_X);
  assert(nmesh_y == NMESH_V);
  io_err = fscanf(input_fp, "%le %le", &(this_run->xmin), &(this_run->xmax));
  io_err = fscanf(input_fp, "%le %le", &(this_run->vmin), &(this_run->vmax));
  for(int32_t ix=0;ix<NMESH_X;ix++) {
    for(int32_t iv=0;iv<NMESH_V;iv++) {
      int32_t ix_, iv_;
      io_err = fscanf(input_fp, "%d %d %le", &ix_, &iv_, &(DF(ix,iv)));
    }
  }
#endif
  
  fclose(input_fp);
}

void output_df_in_run(REAL *df, struct run_param *this_run)
{
  static char filename[256];

  REAL next_output_timing;

  next_output_timing = this_run->output_timing[this_run->output_indx];

  if(next_output_timing < this_run->tnow) {
    sprintf(filename, "%s_%3.1f.dat",this_run->modelname, next_output_timing);
    output_df(df, this_run, filename);
    this_run->output_indx++;
  }
}
