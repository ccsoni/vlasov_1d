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

REAL calc_dtime(struct run_param *this_run, REAL *pot)
{
  REAL dt_x, dt_v;

  dt_x = this_run->delta_x/(fmax(fabs(this_run->vmax), fabs(this_run->vmin)));

  dt_v = FLT_MAX;
#pragma omp parallel for reduction(min:dt_v)
  for(int32_t ix=0;ix<NMESH_X;ix++) {
    REAL acc;

    int32_t ip1 = ix+1; if(ip1>NMESH_X-1) ip1 -= NMESH_X;
    int32_t ip2 = ix+2; if(ip2>NMESH_X-1) ip2 -= NMESH_X;
    int32_t ip3 = ix+3; if(ip3>NMESH_X-1) ip3 -= NMESH_X;
    
    int32_t im1 = ix-1; if(im1<0) im1 += NMESH_X;
    int32_t im2 = ix-2; if(im2<0) im2 += NMESH_X;
    int32_t im3 = ix-3; if(im3<0) im3 += NMESH_X;

#ifdef __POISSON_2ND__
    acc = -0.5*(pot[ip1]-pot[im1])/this_run->delta_x;
#elif __POISSON_4TH__
    acc = ((pot[ip2]-pot[im2])-8.0*(pot[ip1]-pot[im1]))/(12.0*this_run->delta_x);
#elif __POISSON_6TH__
    acc = -((pot[ip3]-pot[im3])-9.0*(pot[ip2]-pot[im2])+45.0*(pot[ip1]-pot[im1]))/(60.0*this_run->delta_x);
#endif

    dt_v = fmin(dt_v, this_run->delta_v/fabs(acc));

  }
  
  REAL dt = 0.1*fmin(dt_x, dt_v);

  REAL next_output_timing = this_run->output_timing[this_run->output_indx];
  if(this_run->tnow + dt > next_output_timing) {
    dt = next_output_timing - this_run->tnow + 1.0e-4;
  }

  this_run->dtime = dt;
  
  return dt;
}

REAL csl5(REAL fuu, REAL fu, REAL f0, REAL fd, REAL fdd,
	  REAL dx, REAL del)
{
  return ((4.0*(fuu-6.5*fu+23.5*f0)+6.0*(9.0*fd-fdd))
               +del*(5.0*(15.0*(f0-fd)-(fu-fdd))
                     +del*((10.0*(3.0*fu-4.0*f0+fd)-5.0*(fuu-fdd))
                           +del*(5.0*(-3.0*(f0-fd)+(fu-fdd))
                                 +del*(6.0*f0-4.0*(fu+fd)+(fuu+fdd))))))/120.0;
}

void calc_flux_x(REAL *df, REAL vel, REAL dtime, struct vflux *flux,
		 struct run_param *this_run, int32_t iv)
{
  REAL del = fabs(vel*dtime/this_run->delta_x);
  REAL lambda = dtime/this_run->delta_x;

  for(int32_t ix=0;ix<NMESH_X;ix++) {
    int32_t ixm3 = ix-3; if(ixm3 < 0)        ixm3 += NMESH_X;
    int32_t ixm2 = ix-2; if(ixm2 < 0)        ixm2 += NMESH_X;
    int32_t ixm1 = ix-1; if(ixm1 < 0)        ixm1 += NMESH_X;
    
    int32_t ixp1 = ix+1; if(ixp1 >= NMESH_X) ixp1 -= NMESH_X;
    int32_t ixp2 = ix+2; if(ixp2 >= NMESH_X) ixp2 -= NMESH_X;
    int32_t ixp3 = ix+3; if(ixp3 >= NMESH_X) ixp3 -= NMESH_X;

    if(vel>0.0) {
      flux[ix].flux_p = csl5(df[ixm2], df[ixm1], df[ix], df[ixp1], df[ixp2],
                             this_run->delta_x, del);
      flux[ix].flux_m = csl5(df[ixm3], df[ixm2], df[ixm1], df[ix], df[ixp1],
                             this_run->delta_x, del);
      flux[ix].flux_p = mp5_bound_L(df[ixm2],df[ixm1],df[ix],df[ixp1],df[ixp2],
                                    flux[ix].flux_p);
      flux[ix].flux_m = mp5_bound_L(df[ixm3],df[ixm2],df[ixm1],df[ix],df[ixp1],
                                    flux[ix].flux_m);
    }else{
      flux[ix].flux_p = csl5(df[ixp3], df[ixp2], df[ixp1], df[ix], df[ixm1],
                             this_run->delta_x, del);
      flux[ix].flux_m = csl5(df[ixp2], df[ixp1], df[ix], df[ixm1], df[ixm2],
                             this_run->delta_x, del);
      flux[ix].flux_p = mp5_bound_R(df[ixm1],df[ix],df[ixp1],df[ixp2],df[ixp3],
                                    flux[ix].flux_p);
      flux[ix].flux_m = mp5_bound_R(df[ixm2],df[ixm1],df[ix],df[ixp1],df[ixp2],
                                    flux[ix].flux_m);
    }    
    flux[ix].flux_p *= vel;
    flux[ix].flux_m *= vel;

  }


#ifdef __PP__
  REAL eps_pp = 0.0;

  for(int32_t ix=0;ix<NMESH_X;ix++) {
    
    int32_t ixp = ix + 1; if(ixp >= NMESH_X) ixp -= NMESH_X;
    int32_t ixm = ix - 1; if(ixm < 0) ixm += NMESH_X;

    REAL theta_p, theta_m, theta;
    REAL dfp, dfm, df_FLp, df_FLm;
    REAL flux_FL;

    /* for flux_p[*] */
    dfp = df[ix]  - 2.0*lambda*flux[ix].flux_p;
    dfm = df[ixp] + 2.0*lambda*flux[ixp].flux_m;

    flux_FL = 0.5*vel*(df[ix]+df[ixp])-0.5*(df[ixp]-df[ix])/lambda;

    theta_p = theta_m = 1.0;
    if(dfp < eps_pp) {
      df_FLp = df[ix] - 2.0*lambda*flux_FL;
      theta_p = (eps_pp - df_FLp)/(dfp - df_FLp + DBL_MIN);
      if(theta_p > 1.0 || theta_p < 0.0) theta_p = 0.0;
    }
    if(dfm < eps_pp) {
      df_FLm = df[ixp] + 2.0*lambda*flux_FL;
      theta_m = (eps_pp - df_FLm)/(dfm - df_FLm + DBL_MIN);
      if(theta_m > 1.0 || theta_m < 0.0) theta_m = 0.0;
    }
    theta = fmin(theta_m, theta_p);
    flux[ix].flux_p = (1.0-theta)*flux_FL + theta*flux[ix].flux_p;
    flux[ixp].flux_m = flux[ix].flux_p;

    /* for flux_m[*] */
    dfp = df[ixm] - 2.0*lambda*flux[ixm].flux_p;
    dfm = df[ix]  + 2.0*lambda*flux[ix].flux_m;

    flux_FL = 0.5*vel*(df[ixm]+df[ix])-0.5*(df[ix]-df[ixm])/lambda;

    theta_p = theta_m = 1.0;
    if(dfp < eps_pp) {
      df_FLp = df[ixm] - 2.0*lambda*flux_FL;
      theta_p = (eps_pp - df_FLp)/(dfp - df_FLp + DBL_MIN);
      if(theta_p > 1.0 || theta_p < 0.0) theta_p = 0.0;      
    }
    if(dfm < eps_pp) {
      df_FLm = df[ix] + 2.0*lambda*flux_FL;
      theta_m = (eps_pp - df_FLm)/(dfm - df_FLm + DBL_MIN);
      if(theta_m > 1.0 || theta_m < 0.0) theta_m = 0.0;      
    }
    theta = fmin(theta_m, theta_p);
    flux[ix].flux_m = (1.0-theta)*flux_FL + theta*flux[ix].flux_m;
    flux[ixm].flux_p = flux[ix].flux_m;

  }
#endif

}

void calc_flux_v(REAL *df, REAL acc, REAL dtime, struct vflux *flux,
		 struct run_param *this_run, int32_t ix)
{
  REAL del = fabs(acc*dtime/this_run->delta_v);
  REAL lambda = dtime/this_run->delta_v;  

  for(int32_t iv=0;iv<NMESH_V;iv++) {
    int32_t ivm3 = iv-3; if(ivm3 < 0)        ivm3 = 0;
    int32_t ivm2 = iv-2; if(ivm2 < 0)        ivm2 = 0;
    int32_t ivm1 = iv-1; if(ivm1 < 0)        ivm1 = 0;
    int32_t ivp1 = iv+1; if(ivp1 >= NMESH_V) ivp1 = NMESH_V-1;
    int32_t ivp2 = iv+2; if(ivp2 >= NMESH_V) ivp2 = NMESH_V-1;
    int32_t ivp3 = iv+3; if(ivp3 >= NMESH_V) ivp3 = NMESH_V-1;

    if(acc > 0.0) {
      flux[iv].flux_p = csl5(df[ivm2], df[ivm1], df[iv], df[ivp1], df[ivp2],
                             this_run->delta_v, del);
      flux[iv].flux_m = csl5(df[ivm3], df[ivm2], df[ivm1], df[iv], df[ivp1],
                             this_run->delta_v, del);
      flux[iv].flux_p = mp5_bound_L(df[ivm2],df[ivm1],df[iv],df[ivp1],df[ivp2],
                                    flux[iv].flux_p);
      flux[iv].flux_m = mp5_bound_L(df[ivm3],df[ivm2],df[ivm1],df[iv],df[ivp1],
                                    flux[iv].flux_m);
    }else{
      flux[iv].flux_p = csl5(df[ivp3], df[ivp2], df[ivp1], df[iv], df[ivm1],
                             this_run->delta_v, del);
      flux[iv].flux_m = csl5(df[ivp2], df[ivp1], df[iv], df[ivm1], df[ivm2],
                             this_run->delta_v, del);
      flux[iv].flux_p = mp5_bound_R(df[ivm1],df[iv],df[ivp1],df[ivp2],df[ivp3],
                                    flux[iv].flux_p);
      flux[iv].flux_m = mp5_bound_R(df[ivm2],df[ivm1],df[iv],df[ivp1],df[ivp2],
                                    flux[iv].flux_m);
    }
    flux[iv].flux_p *= acc;
    flux[iv].flux_m *= acc;

  }

#ifdef __PP__
  REAL eps_pp = 0.0;

  for(int32_t iv=0;iv<NMESH_V;iv++) {
    
    int32_t ivp = iv + 1; if(ivp >= NMESH_V) ivp = NMESH_V-1;
    int32_t ivm = iv - 1; if(ivm < 0) ivm = 0;

    REAL theta_p, theta_m, theta;
    REAL dfp, dfm, df_FLp, df_FLm;
    REAL flux_FL;

    /* for flux_p[*] */
    dfp = df[iv]  - 2.0*lambda*flux[iv].flux_p;
    dfm = df[ivp] + 2.0*lambda*flux[ivp].flux_m;

    flux_FL = 0.5*acc*(df[iv]+df[ivp])-0.5*(df[ivp]-df[iv])/lambda;

    theta_p = theta_m = 1.0;
    if(dfp < eps_pp) {
      df_FLp = df[iv] - 2.0*lambda*flux_FL;
      theta_p = (eps_pp - df_FLp)/(dfp - df_FLp + DBL_MIN);
      if(theta_p > 1.0 || theta_p < 0.0) theta_p = 0.0;
    }
    if(dfm < eps_pp) {
      df_FLm = df[ivp] + 2.0*lambda*flux_FL;
      theta_m = (eps_pp - df_FLm)/(dfm - df_FLm + DBL_MIN);
      if(theta_m > 1.0 || theta_m < 0.0) theta_m = 0.0;
    }
    theta = fmin(theta_m, theta_p);
    flux[iv].flux_p = (1.0-theta)*flux_FL + theta*flux[iv].flux_p;
    flux[ivp].flux_m = flux[iv].flux_p;

    /* for flux_m[*] */
    dfp = df[ivm] - 2.0*lambda*flux[ivm].flux_p;
    dfm = df[iv]  + 2.0*lambda*flux[iv].flux_m;

    flux_FL = 0.5*acc*(df[ivm]+df[iv])-0.5*(df[iv]-df[ivm])/lambda;

    theta_p = theta_m = 1.0;
    if(dfp < eps_pp) {
      df_FLp = df[ivm] - 2.0*lambda*flux_FL;
      theta_p = (eps_pp - df_FLp)/(dfp - df_FLp + DBL_MIN);
      if(theta_p > 1.0 || theta_p < 0.0) theta_p = 0.0;
    }
    if(dfm < eps_pp) {
      df_FLm = df[iv] + 2.0*lambda*flux_FL;
      theta_m = (eps_pp - df_FLm)/(dfm - df_FLm + DBL_MIN);
      if(theta_m > 1.0 || theta_m < 0.0) theta_m = 0.0;
    }
    theta = fmin(theta_m, theta_p);
    flux[iv].flux_m = (1.0-theta)*flux_FL + theta*flux[iv].flux_m;
    flux[ivm].flux_p = flux[iv].flux_m;
  }
#endif // __PP__

}

void integrate_x(REAL *df, REAL dtime, struct run_param *this_run)
{
  REAL lambda = dtime/this_run->delta_x;
  
#pragma omp parallel for schedule(auto)
  for(int32_t iv=0;iv<NMESH_V;iv++) {
    struct vflux *flux = (struct vflux *)malloc(sizeof(struct vflux)*NMESH_X);
    REAL *df_copy = (REAL *)malloc(sizeof(REAL)*NMESH_X);

    REAL vel = this_run->vmin + this_run->delta_v*((REAL)iv+0.5);

    for(int32_t ix=0;ix<NMESH_X;ix++) df_copy[ix] = DF(ix,iv);

    calc_flux_x(df_copy, vel, dtime, flux, this_run,iv);

    for(int32_t ix=0;ix<NMESH_X;ix++) {
      df_copy[ix] -= lambda*(flux[ix].flux_p - flux[ix].flux_m);
    }

    for(int32_t ix=0;ix<NMESH_X;ix++) DF(ix,iv) = df_copy[ix];

    free(flux);
    free(df_copy);
  }
}

void integrate_v(REAL *df, REAL *pot, REAL dtime, struct run_param *this_run)
{
  REAL lambda = dtime/this_run->delta_x;
  
#pragma omp parallel for schedule(auto)
  for(int32_t ix=0;ix<NMESH_X;ix++) {
    struct vflux *flux = (struct vflux *)malloc(sizeof(struct vflux)*NMESH_V);
    REAL *df_copy = (REAL *)malloc(sizeof(REAL)*NMESH_V);

    int32_t ixp3 = ix + 3; if(ixp3 >= NMESH_X) ixp3 -= NMESH_X;
    int32_t ixp2 = ix + 2; if(ixp2 >= NMESH_X) ixp2 -= NMESH_X;
    int32_t ixp1 = ix + 1; if(ixp1 >= NMESH_X) ixp1 -= NMESH_X;

    int32_t ixm1 = ix - 1; if(ixm1 < 0) ixm1 += NMESH_X;
    int32_t ixm2 = ix - 2; if(ixm2 < 0) ixm2 += NMESH_X;
    int32_t ixm3 = ix - 3; if(ixm3 < 0) ixm3 += NMESH_X;

#ifdef __POISSON_2ND__
    REAL acc = -0.5*(pot[ixp1] - pot[ixm1])/this_run->delta_x;
#elif __POISSON_4TH__
    REAL acc = ((pot[ixp2]-pot[ixm2])-8.0*(pot[ixp1]-pot[ixm1]))/(12.0*this_run->delta_x);
#elif __POISSON_6TH__
    REAL acc = -((pot[ixp3]-pot[ixm3])-9.0*(pot[ixp2]-pot[ixm2])+45.0*(pot[ixp1]-pot[ixm1]))/(60.0*this_run->delta_x);
#endif

    for(int iv=0;iv<NMESH_V;iv++) df_copy[iv] = DF(ix,iv);

    calc_flux_v(df_copy, acc, dtime, flux, this_run, ix);
    for(int iv=0;iv<NMESH_V;iv++) {
      df_copy[iv] -= lambda*(flux[iv].flux_p-flux[iv].flux_m);
    }

    for(int iv=0;iv<NMESH_V;iv++) DF(ix,iv) = df_copy[iv];
    
    free(df_copy);
    free(flux);
  }
}
