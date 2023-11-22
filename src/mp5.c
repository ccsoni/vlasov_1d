#include <math.h>

#include "run_param.h"

REAL minmod2(REAL x, REAL y)
{
  REAL minmod;
  REAL abs_val;
  REAL min_val;
  
  min_val = fmin(fabs(x),fabs(y));
  minmod = 0.5*(copysign(1.0,x)+copysign(1.0,y))*min_val;

  return minmod;
}

REAL minmod4(REAL w, REAL x, REAL y, REAL z)
{
  REAL minmod;
  REAL abs_val;
  REAL min_val;

  min_val = fmin(fabs(w),fmin(fabs(x),fmin(fabs(y),fabs(z))));
  abs_val = (copysign(1.0,w)+copysign(1.0,y))*(copysign(1.0,w)+copysign(1.0,z));
  minmod = 0.125*(copysign(1.0,w)+copysign(1.0,x))*fabs(abs_val)*min_val;

  return minmod;
}

REAL mp5_R(REAL fjm2, REAL fjm1, REAL fj, REAL fjp1, REAL fjp2)
{
  REAL f_int, f_mp, f_R;
  REAL f_ul, f_av, f_md, f_lc, f_min, f_max;
  REAL C1 = 1.0/60.0;
  REAL C2 = 4.0/3.0;

  REAL ddfjm1, ddfj, ddfjp1;
  REAL ddfjph, ddfjmh;

#define EPS (0.0)
#define ALPHA (4.0)

  f_int = C1*(2.0*fjp2 - 13.0*fjp1 + 47.0*fj + 27.0*fjm1 - 3.0*fjm2);
  f_mp = fj + minmod2(fjm1-fj, ALPHA*(fj-fjp1));
  if((f_int-fj)*(f_int-f_mp) < EPS) {
    f_R = f_int;
  }else{
    ddfjp1 = fjp2 - 2.0*fjp1 + fj;
    ddfj   = fjp1 - 2.0*fj   + fjm1;
    ddfjm1 = fj   - 2.0*fjm1 + fjm2;
    ddfjph  = minmod4(4.0*ddfj-ddfjp1, 4.0*ddfjp1-ddfj, ddfj, ddfjp1);
    ddfjmh  = minmod4(4.0*ddfj-ddfjm1, 4.0*ddfjm1-ddfj, ddfj, ddfjm1);
    
    f_ul = fj + ALPHA*(fj-fjp1);
    f_av = 0.5*(fj+fjm1);
    f_md = f_av - 0.5*ddfjmh;

    f_lc = fj + 0.5*(fj-fjp1) + C2*ddfjph;
#if 0
    f_ul = MAX(f_ul, 0.0);
    f_md = MAX(f_md, 0.0);
    f_lc = MAX(f_lc, 0.0);
#endif

    f_min = fmax(fmin(fj,fmin(fjm1,f_md)), 
		 fmin(fj,fmin(f_ul,f_lc)));
    f_max = fmin(fmax(fj,fmax(fjm1,f_md)),
		 fmax(fj,fmax(f_ul,f_lc)));
    
    f_R =  f_int + minmod2(f_min-f_int, f_max-f_int);

  }

#undef ALPHA
#undef EPS

  return f_R;

}

REAL mp5_L(REAL fjm2, REAL fjm1, REAL fj, REAL fjp1, REAL fjp2)
{
  REAL f_int, f_mp, f_L;
  REAL f_ul, f_av, f_md, f_lc, f_min, f_max;
  REAL C1 = 1.0/60.0;
  REAL C2 = 4.0/3.0;

  REAL ddfjm1, ddfj, ddfjp1;
  REAL ddfjph, ddfjmh;

#define EPS (0.0)
#define ALPHA (4.0)

  f_int = C1*(2.0*fjm2 - 13.0*fjm1 + 47.0*fj + 27.0*fjp1 - 3.0*fjp2);
  f_mp = fj + minmod2(fjp1-fj, ALPHA*(fj-fjm1));
  if((f_int-fj)*(f_int-f_mp) < EPS) {
    f_L = f_int;
  }else{
    ddfjm1 = fjm2 - 2.0*fjm1 + fj;
    ddfj   = fjm1 - 2.0*fj   + fjp1;
    ddfjp1 = fj   - 2.0*fjp1 + fjp2;
    ddfjph  = minmod4(4.0*ddfj-ddfjp1, 4.0*ddfjp1-ddfj, ddfj, ddfjp1);
    ddfjmh  = minmod4(4.0*ddfj-ddfjm1, 4.0*ddfjm1-ddfj, ddfj, ddfjm1);
    
    f_ul = fj + ALPHA*(fj-fjm1);
    f_av = 0.5*(fj+fjp1);
    f_md = f_av - 0.5*ddfjph;

    f_lc = fj + 0.5*(fj-fjm1) + C2*ddfjmh;
#if 0
    f_ul = MAX(f_ul, 0.0);
    f_md = MAX(f_md, 0.0);
    f_lc = MAX(f_lc, 0.0);
#endif
    
    f_min = fmax(fmin(fj,fmin(fjp1,f_md)), 
		 fmin(fj,fmin(f_ul,f_lc)));
    f_max = fmin(fmax(fj,fmax(fjp1,f_md)),
		 fmax(fj,fmax(f_ul,f_lc)));
    
    f_L =  f_int + minmod2(f_min-f_int, f_max-f_int);
  }

#undef ALPHA
#undef EPS

  return f_L;
}

REAL mp5_bound_L(REAL fjm2, REAL fjm1, REAL fj, REAL fjp1, REAL fjp2, REAL f_in)
{
  REAL f_mp, f_L;
  REAL f_ul, f_av, f_md, f_lc, f_min, f_max;
  REAL C1 = 1.0/60.0;
  REAL C2 = 4.0/3.0;

  REAL ddfjm1, ddfj, ddfjp1;
  REAL ddfjph, ddfjmh;

#define EPS (0.0)
#define ALPHA (4.0)

  f_mp = fj + minmod2(fjp1-fj, ALPHA*(fj-fjm1));
  if((f_in-fj)*(f_in-f_mp) < EPS) {
    f_L = f_in;
  }else{
    ddfjm1 = fjm2 - 2.0*fjm1 + fj;
    ddfj   = fjm1 - 2.0*fj   + fjp1;
    ddfjp1 = fj   - 2.0*fjp1 + fjp2;
    ddfjph  = minmod4(4.0*ddfj-ddfjp1, 4.0*ddfjp1-ddfj, ddfj, ddfjp1);
    ddfjmh  = minmod4(4.0*ddfj-ddfjm1, 4.0*ddfjm1-ddfj, ddfj, ddfjm1);
    
    f_ul = fj + ALPHA*(fj-fjm1);
    f_av = 0.5*(fj+fjp1);
    f_md = f_av - 0.5*ddfjph;

    f_lc = fj + 0.5*(fj-fjm1) + C2*ddfjmh;
#if 0
    f_ul = MAX(f_ul, 0.0);
    f_md = MAX(f_md, 0.0);
    f_lc = MAX(f_lc, 0.0);
#endif
    
    f_min = fmax(fmin(fj,fmin(fjp1,f_md)), 
		 fmin(fj,fmin(f_ul,f_lc)));
    f_max = fmin(fmax(fj,fmax(fjp1,f_md)),
		 fmax(fj,fmax(f_ul,f_lc)));
    
    f_L =  f_in + minmod2(f_min-f_in, f_max-f_in);
  }  

#undef EPS
#undef ALPHA  

  return f_L;
}

REAL mp5_bound_R(REAL fjm2, REAL fjm1, REAL fj, REAL fjp1, REAL fjp2, REAL f_in)
{
  REAL f_mp, f_R;
  REAL f_ul, f_av, f_md, f_lc, f_min, f_max;
  REAL C1 = 1.0/60.0;
  REAL C2 = 4.0/3.0;

  REAL ddfjm1, ddfj, ddfjp1;
  REAL ddfjph, ddfjmh;

#define EPS (0.0)
#define ALPHA (4.0)

  f_mp = fj + minmod2(fjm1-fj, ALPHA*(fj-fjp1));
  if((f_in-fj)*(f_in-f_mp) < EPS) {
    f_R = f_in;
  }else{
    ddfjp1 = fjp2 - 2.0*fjp1 + fj;
    ddfj   = fjp1 - 2.0*fj   + fjm1;
    ddfjm1 = fj   - 2.0*fjm1 + fjm2;
    ddfjph  = minmod4(4.0*ddfj-ddfjp1, 4.0*ddfjp1-ddfj, ddfj, ddfjp1);
    ddfjmh  = minmod4(4.0*ddfj-ddfjm1, 4.0*ddfjm1-ddfj, ddfj, ddfjm1);
    
    f_ul = fj + ALPHA*(fj-fjp1);
    f_av = 0.5*(fj+fjm1);
    f_md = f_av - 0.5*ddfjmh;

    f_lc = fj + 0.5*(fj-fjp1) + C2*ddfjph;
#if 0
    f_ul = MAX(f_ul, 0.0);
    f_md = MAX(f_md, 0.0);
    f_lc = MAX(f_lc, 0.0);
#endif

    f_min = fmax(fmin(fj,fmin(fjm1,f_md)), 
		 fmin(fj,fmin(f_ul,f_lc)));
    f_max = fmin(fmax(fj,fmax(fjm1,f_md)),
		 fmax(fj,fmax(f_ul,f_lc)));
    
    f_R =  f_in + minmod2(f_min-f_in, f_max-f_in);

  }

#undef ALPHA
#undef EPS

  return f_R;

}
