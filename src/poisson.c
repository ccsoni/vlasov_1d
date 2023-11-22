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

void calc_dens_field(REAL *df, REAL *dens, struct run_param *this_run)
{
#pragma omp parallel for
  for(int32_t ix=0;ix<NMESH_X;ix++) {
    dens[ix] = 0.0;
    for(int32_t iv=0;iv<NMESH_V;iv++) {
      dens[ix] += DF(ix,iv)*this_run->delta_v;
    }
  }
}

void calc_poten(REAL *dens, REAL *pot, struct run_param *this_run)
{
  static REAL work_real[NMESH_X], work_imag[NMESH_X];
  static REAL data_real[NMESH_X], data_imag[NMESH_X];
  static REAL green[NGREEN];

  REAL dx = this_run->delta_x;

#pragma omp parallel for
  for(int ix=0;ix<NMESH_X;ix++) dens[ix] *= dx;

#pragma omp parallel for
#ifdef __POISSON_2ND__
  for(int k=1;k<NGREEN;k++) {
    REAL sinc;
    sinc = sin(M_PI*(REAL)k*dx)/dx;
    green[k]=-M_PI/(sinc*sinc);
  }
  green[0] = 0.0;
#elif __POISSON_4TH__
  for(int k=1;k<NGREEN;k++) {
    green[k]=12.0*M_PI*SQR(dx)/(SQR(sin(2.0*M_PI*(REAL)k*dx))-16.0*SQR(sin(M_PI*(REAL)k*dx)));
  }
  green[0] = 0.0;
#elif __POISSON_6TH__
  for(int k=1;k<NGREEN;k++) {
    green[k]=-180.0*M_PI*SQR(dx)/(2.0*SQR(sin(3.0*M_PI*(REAL)k*dx))-27.0*SQR(sin(2.0*M_PI*(REAL)k*dx))+270.0*SQR(sin(M_PI*(REAL)k*dx)));
  }
  green[0] = 0.0;
#endif


#pragma omp parallel for
  for(int i=0;i<NMESH_X;i++) {
    work_real[i]=0.0;
    work_imag[i]=0.0;
    for(int j=0;j<NMESH_X;j++) {
      work_real[i] += (dens[j]*cos(2.0*M_PI*(REAL)(i*j)/(REAL)(NMESH_X)));
      work_imag[i] -= (dens[j]*sin(2.0*M_PI*(REAL)(i*j)/(REAL)(NMESH_X)));
    }
  }

#pragma omp parallel for
  for(int k=0;k<NMESH_X/2;k++) {
    work_real[k] *= green[k];
    work_imag[k] *= green[k];
    work_real[NMESH_X/2+k] *= green[NMESH_X/2-k];
    work_imag[NMESH_X/2+k] *= green[NMESH_X/2-k];
  }

#pragma omp parallel for
  for(int i=0;i<NMESH_X;i++) {
    data_real[i] = 0.0;
    data_imag[i] = 0.0;
    for(int j=0;j<NMESH_X;j++) {
      data_real[i] += (work_real[j]*cos(2.0*M_PI*(REAL)(i*j)/(REAL)(NMESH_X))-
                       work_imag[j]*sin(2.0*M_PI*(REAL)(i*j)/(REAL)(NMESH_X)));
      data_imag[i] += (work_real[j]*sin(2.0*M_PI*(REAL)(i*j)/(REAL)(NMESH_X))+
                       work_imag[j]*cos(2.0*M_PI*(REAL)(i*j)/(REAL)(NMESH_X)));
    }
  }

#pragma omp parallel for
  for(int i=0;i<NMESH_X;i++) pot[i] = data_real[i];
  
}
