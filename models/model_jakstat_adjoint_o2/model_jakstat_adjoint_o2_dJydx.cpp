
#include <include/symbolic_functions.h>
#include <string.h>
#include <include/udata.h>
#include <include/udata_accessors.h>
#include "model_jakstat_adjoint_o2_w.h"

int dJydx_model_jakstat_adjoint_o2(realtype t, int it, realtype *dJydx, realtype *y, N_Vector x, realtype *dydx, realtype *my, realtype *sigma_y, void *user_data) {
int status = 0;
UserData *udata = (UserData*) user_data;
realtype *x_tmp = N_VGetArrayPointer(x);
status = w_model_jakstat_adjoint_o2(t,x,NULL,user_data);
int iy;
if(!amiIsNaN(my[0*nt+it])){
    iy = 0;
  dJydx[it+(1+0*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  dJydx[it+(1+1*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*3]*1.0;
  dJydx[it+(1+2*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*6]*1.0;
  dJydx[it+(1+3*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*9]*1.0;
  dJydx[it+(1+4*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*12]*1.0;
  dJydx[it+(1+5*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*15]*1.0-dydx[69]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*5.0E-1;
  dJydx[it+(1+6*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*18]*1.0;
  dJydx[it+(1+7*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*21]*1.0;
  dJydx[it+(1+8*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*24]*1.0;
  dJydx[it+(1+9*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*27]*1.0;
  dJydx[it+(1+10*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*30]*1.0;
  dJydx[it+(1+11*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*33]*1.0;
  dJydx[it+(1+12*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*36]*1.0;
  dJydx[it+(1+13*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*39]*1.0;
  dJydx[it+(1+14*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*42]*1.0-dydx[96]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*5.0E-1;
  dJydx[it+(1+15*9)*nt] += dydx[54]*(1.0/(sigma_y[0]*sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*1.0+1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*45]*1.0);
  dJydx[it+(1+16*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*48]*1.0;
  dJydx[it+(1+17*9)*nt] += dydx[54]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*51]*1.0;
  dJydx[it+(2+0*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*-5.0E-1;
  dJydx[it+(2+1*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*3]*1.0;
  dJydx[it+(2+2*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*6]*1.0;
  dJydx[it+(2+3*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*9]*1.0;
  dJydx[it+(2+4*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*12]*1.0;
  dJydx[it+(2+5*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*15]*1.0-dydx[123]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*5.0E-1;
  dJydx[it+(2+6*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*18]*1.0;
  dJydx[it+(2+7*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*21]*1.0;
  dJydx[it+(2+8*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*24]*1.0;
  dJydx[it+(2+9*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*27]*1.0;
  dJydx[it+(2+10*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*30]*1.0;
  dJydx[it+(2+11*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*33]*1.0;
  dJydx[it+(2+12*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*36]*1.0;
  dJydx[it+(2+13*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*39]*1.0;
  dJydx[it+(2+14*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*42]*1.0-dydx[150]*1.0/(sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*5.0E-1;
  dJydx[it+(2+15*9)*nt] += dydx[108]*(1.0/(sigma_y[0]*sigma_y[0]*sigma_y[0])*(my[it+nt*0]*2.0-y[it+nt*0]*2.0)*1.0+1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*45]*1.0);
  dJydx[it+(2+16*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*48]*1.0;
  dJydx[it+(2+17*9)*nt] += dydx[108]*1.0/(sigma_y[0]*sigma_y[0])*y[it+nt*51]*1.0;
}
if(!amiIsNaN(my[1*nt+it])){
    iy = 1;
  dJydx[it+(0+0*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  dJydx[it+(0+1*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*4]*1.0;
  dJydx[it+(0+2*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*7]*1.0;
  dJydx[it+(0+3*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*10]*1.0;
  dJydx[it+(0+4*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*13]*1.0;
  dJydx[it+(0+5*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*16]*1.0-dydx[16]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*5.0E-1;
  dJydx[it+(0+6*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*19]*1.0;
  dJydx[it+(0+7*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*22]*1.0;
  dJydx[it+(0+8*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*25]*1.0;
  dJydx[it+(0+9*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*28]*1.0;
  dJydx[it+(0+10*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*31]*1.0;
  dJydx[it+(0+11*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*34]*1.0;
  dJydx[it+(0+12*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*37]*1.0;
  dJydx[it+(0+13*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*40]*1.0-dydx[40]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*5.0E-1;
  dJydx[it+(0+14*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*43]*1.0;
  dJydx[it+(0+15*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*46]*1.0;
  dJydx[it+(0+16*9)*nt] += dydx[1]*(1.0/(sigma_y[1]*sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*1.0+1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*49]*1.0);
  dJydx[it+(0+17*9)*nt] += dydx[1]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*52]*1.0;
  dJydx[it+(1+0*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  dJydx[it+(1+1*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*4]*1.0;
  dJydx[it+(1+2*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*7]*1.0;
  dJydx[it+(1+3*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*10]*1.0;
  dJydx[it+(1+4*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*13]*1.0;
  dJydx[it+(1+5*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*16]*1.0-dydx[70]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*5.0E-1;
  dJydx[it+(1+6*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*19]*1.0;
  dJydx[it+(1+7*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*22]*1.0;
  dJydx[it+(1+8*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*25]*1.0;
  dJydx[it+(1+9*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*28]*1.0;
  dJydx[it+(1+10*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*31]*1.0;
  dJydx[it+(1+11*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*34]*1.0;
  dJydx[it+(1+12*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*37]*1.0;
  dJydx[it+(1+13*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*40]*1.0-dydx[94]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*5.0E-1;
  dJydx[it+(1+14*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*43]*1.0;
  dJydx[it+(1+15*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*46]*1.0;
  dJydx[it+(1+16*9)*nt] += dydx[55]*(1.0/(sigma_y[1]*sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*1.0+1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*49]*1.0);
  dJydx[it+(1+17*9)*nt] += dydx[55]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*52]*1.0;
  dJydx[it+(2+0*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*-5.0E-1;
  dJydx[it+(2+1*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*4]*1.0;
  dJydx[it+(2+2*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*7]*1.0;
  dJydx[it+(2+3*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*10]*1.0;
  dJydx[it+(2+4*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*13]*1.0;
  dJydx[it+(2+5*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*16]*1.0-dydx[124]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*5.0E-1;
  dJydx[it+(2+6*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*19]*1.0;
  dJydx[it+(2+7*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*22]*1.0;
  dJydx[it+(2+8*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*25]*1.0;
  dJydx[it+(2+9*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*28]*1.0;
  dJydx[it+(2+10*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*31]*1.0;
  dJydx[it+(2+11*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*34]*1.0;
  dJydx[it+(2+12*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*37]*1.0;
  dJydx[it+(2+13*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*40]*1.0-dydx[148]*1.0/(sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*5.0E-1;
  dJydx[it+(2+14*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*43]*1.0;
  dJydx[it+(2+15*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*46]*1.0;
  dJydx[it+(2+16*9)*nt] += dydx[109]*(1.0/(sigma_y[1]*sigma_y[1]*sigma_y[1])*(my[it+nt*1]*2.0-y[it+nt*1]*2.0)*1.0+1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*49]*1.0);
  dJydx[it+(2+17*9)*nt] += dydx[109]*1.0/(sigma_y[1]*sigma_y[1])*y[it+nt*52]*1.0;
}
if(!amiIsNaN(my[2*nt+it])){
    iy = 2;
}
return(status);

}


