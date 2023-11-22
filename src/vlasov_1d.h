#pragma once

#include <stdbool.h>
#include "run_param.h"

struct vlasov_flux {
  REAL flux_xp, flux_xm;
  REAL flux_vp, flux_vm;
};

struct vflux {
  REAL flux_p, flux_m;
  bool pp_correction;
};
