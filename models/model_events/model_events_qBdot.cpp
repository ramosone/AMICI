
#include "amici/symbolic_functions.h"
#include "amici/defines.h" //realtype definition
typedef amici::realtype realtype;
#include <cmath> 

using namespace amici;

void qBdot_model_events(realtype *qBdot, const int ip, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdp) {
switch (ip) {
  case 0: {
  qBdot[0] = -x[0]*xB[0]*(h[3]-1.0);

  } break;

  case 1: {
  qBdot[0] = -x[0]*xB[1]*exp(t*(-1.0/1.0E1));

  } break;

  case 2: {
  qBdot[0] = x[1]*xB[1];

  } break;

}
}

