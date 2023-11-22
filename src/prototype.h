#pragma once

#include "run_param.h"

void setup_IC(REAL*, struct run_param*);
void input_df(REAL*, struct run_param*, char*);
void output_df(REAL*, struct run_param*, char*);
void output_df_in_run(REAL*, struct run_param*);
void output_diag(REAL*, struct run_param*, REAL*, REAL*);

REAL mp5_bound_R(REAL, REAL, REAL, REAL, REAL, REAL);
REAL mp5_bound_L(REAL, REAL, REAL, REAL, REAL, REAL);
void integrate_x(REAL*, REAL, struct run_param*);
void integrate_v(REAL*, REAL*, REAL, struct run_param*);
void calc_dens_field(REAL*, REAL*, struct run_param*);
void calc_poten(REAL*, REAL*, struct run_param*);
REAL calc_dtime(struct run_param*, REAL *);

